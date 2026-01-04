package main

import (
	"fmt"
	"math"
)

type Box struct{ W, H, L int } // mm

type ItemType struct {
	Name  string // item type name (identifier for printing/debugging)
	Count int    // how many items of this type remain to be packed
	Dim   Box    // original item dimensions (in mm)
}

type Container struct {
	Name string // container type name (identifier for printing/debugging)
	Dim  Box    // container inner dimensions (in mm)
}

// Orientation is an axis-aligned permutation of an item's dimensions.
// It represents how the item is rotated relative to container axes.
type Orientation struct{ W, H, L int } // aligned to container axes

type Placement struct {
	InBox   Box         // the sub-box (region) we packed into (size only; no position tracking)
	Item    string      // item type name that was placed
	Ori     Orientation // chosen orientation for this placement
	Count   int         // how many items were placed in this placement
	Nx, Ny  int         // grid count along W (Nx) and H (Ny) in the cross-section
	NLayers int         // number of layers along L used by this placement
	Used    Box         // used block size inside InBox (Nx*W, Ny*H, NLayers*L)
}

func cmToMm(x float64) int { return int(math.Round(x * 10)) } // cm -> mm (rounded)

func perm3(a, b, c int) [][3]int {
	return [][3]int{
		{a, b, c}, {a, c, b},
		{b, a, c}, {b, c, a},
		{c, a, b}, {c, b, a},
	}
}

// AllOrientations returns all unique orientations of the item dimensions,
// aligned to container axes (W, H, L).
func AllOrientations(it ItemType) []Orientation {
	seen := map[[3]int]struct{}{}    // set of (W,H,L) to remove duplicates (e.g., cubes)
	out := make([]Orientation, 0, 6) // up to 6 unique permutations

	for _, p := range perm3(it.Dim.W, it.Dim.H, it.Dim.L) { // generate all 3D permutations
		o := Orientation{W: p[0], H: p[1], L: p[2]} // treat permuted dims as axis-aligned orientation
		key := [3]int{o.W, o.H, o.L}                // de-duplication key
		if _, ok := seen[key]; ok {
			continue
		}
		seen[key] = struct{}{} // mark orientation as seen
		out = append(out, o)   // collect unique orientation
	}

	return out // all unique axis-aligned orientations
}

func vol(b Box) int64 { return int64(b.W) * int64(b.H) * int64(b.L) } // volume in mm^3

func ceilDiv(a, b int) int {
	if b <= 0 {
		return 0
	}
	return (a + b - 1) / b
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

type bestCandidate struct {
	itemIdx     int         // index of chosen item type in items slice
	ori         Orientation // chosen orientation for this candidate
	countPlaced int         // how many items we can place (bounded by Count and capacity)
	nx, ny      int         // grid fit along W/H
	nLayers     int         // layers used along L
	used        Box         // block size used by this candidate (in mm)
	unusedVol   int64       // heuristic: (box volume - placed*item volume)
}

func betterCandidate(a, b bestCandidate) bool {
	if a.countPlaced != b.countPlaced {
		return a.countPlaced > b.countPlaced
	}
	if a.unusedVol != b.unusedVol {
		return a.unusedVol < b.unusedVol
	}
	// tie-breaker: prefer larger used volume (more compact leftovers)
	ua := vol(a.used)
	ub := vol(b.used)
	if ua != ub {
		return ua > ub
	}
	return a.nLayers > b.nLayers
}

func copyItems(items []ItemType) []ItemType {
	out := make([]ItemType, len(items)) // allocate a copy to avoid mutating input
	copy(out, items)                    // copy item records (including Count)
	return out                          // independent slice copy
}

func allPlaced(items []ItemType) bool {
	for _, it := range items { // scan all item types
		if it.Count > 0 { // any remaining items means "not done"
			return false
		}
	}
	return true
}

func sumPlaced(ps []Placement) int {
	total := 0
	for _, p := range ps { // sum per-placement counts
		total += p.Count // accumulate placed items
	}
	return total
}

func containerVol(c Container) int64 { return vol(c.Dim) } // container volume in mm^3 (used for tie-breaking)

// packBoxRecursive greedily packs items into 'box' using axis-aligned grid placements,
// then recursively tries to pack remaining items into leftover sub-boxes.
//
// Leftover partition is a simple non-overlapping guillotine split:
// - right: leftover along W for the used length region
// - top: leftover along H for the used length region (only over usedW)
// - back: leftover along L for the full cross-section
func packBoxRecursive(box Box, items []ItemType) ([]Placement, []ItemType) {
	// Find best single placement in this box.
	found := false         // whether any candidate placement was found for this box
	var best bestCandidate // best candidate found so far

	for i := range items { // iterate over item types
		if items[i].Count <= 0 { // skip items with nothing left to pack
			continue
		}
		for _, o := range AllOrientations(items[i]) { // try all axis-aligned orientations
			if o.W <= 0 || o.H <= 0 || o.L <= 0 { // guard against invalid dimensions
				continue
			}
			if o.W > box.W || o.H > box.H || o.L > box.L { // must fit within current box
				continue
			}

			nx := box.W / o.W                  // how many items fit along width
			ny := box.H / o.H                  // how many items fit along height
			nz := box.L / o.L                  // how many items fit along length (max layers)
			if nx <= 0 || ny <= 0 || nz <= 0 { // must fit at least 1 in each dimension
				continue
			}

			capacity := nx * ny * nz                   // maximum items for this orientation in this box
			placed := minInt(items[i].Count, capacity) // actual items placed limited by remaining Count
			if placed <= 0 {                           // nothing to place
				continue
			}

			perLayer := nx * ny                              // items per single layer along length
			nLayers := minInt(nz, ceilDiv(placed, perLayer)) // layers actually needed to place 'placed' items
			if nLayers <= 0 {                                // should not happen if placed>0, but keep safe
				continue
			}

			used := Box{W: nx * o.W, H: ny * o.H, L: nLayers * o.L}                // bounding block we occupy
			unusedVol := vol(box) - int64(placed)*vol(Box{W: o.W, H: o.H, L: o.L}) // heuristic leftover volume
			if unusedVol < 0 {
				unusedVol = 0
			}

			cand := bestCandidate{ // candidate summary used for comparing options
				itemIdx:     i,
				ori:         o,
				countPlaced: placed,
				nx:          nx,
				ny:          ny,
				nLayers:     nLayers,
				used:        used,
				unusedVol:   unusedVol,
			}

			if !found || betterCandidate(cand, best) { // choose the best candidate by our comparator
				found = true
				best = cand
			}
		}
	}

	if !found {
		return nil, items // nothing fits into this box
	}

	// Apply placement: decrement count.
	items[best.itemIdx].Count -= best.countPlaced // update remaining Count for the chosen item type
	if items[best.itemIdx].Count < 0 {
		items[best.itemIdx].Count = 0
	}

	placements := []Placement{{
		InBox:   box,                      // box we packed into
		Item:    items[best.itemIdx].Name, // item type name
		Ori:     best.ori,                 // orientation chosen
		Count:   best.countPlaced,         // how many were placed
		Nx:      best.nx,                  // grid along W
		Ny:      best.ny,                  // grid along H
		NLayers: best.nLayers,             // layers along L
		Used:    best.used,                // occupied block dimensions
	}}

	// Build leftover sub-boxes (non-overlapping).
	usedW, usedH, usedL := best.used.W, best.used.H, best.used.L // occupied extents for splitting
	right := Box{W: box.W - usedW, H: box.H, L: usedL}           // leftover region to the "right" (along W)
	top := Box{W: usedW, H: box.H - usedH, L: usedL}             // leftover region on "top" (along H)
	back := Box{W: box.W, H: box.H, L: box.L - usedL}            // leftover region "behind" (along L)

	sub := make([]Box, 0, 3) // list of leftover sub-boxes we will recurse into
	if right.W > 0 && right.H > 0 && right.L > 0 {
		sub = append(sub, right) // keep only positive-volume leftovers
	}
	if top.W > 0 && top.H > 0 && top.L > 0 {
		sub = append(sub, top)
	}
	if back.W > 0 && back.H > 0 && back.L > 0 {
		sub = append(sub, back)
	}

	// If the last layer is only partially filled, expose the internal empty cells of that last layer
	// as additional leftover sub-boxes (so we can "finish" the layer with other items).
	//
	// IMPORTANT: We do not track coordinates, only sizes. These sub-boxes are still non-overlapping
	// with right/top/back because they are inside the used cross-section (usedW×usedH) and only
	// within the last layer thickness (best.ori.L).
	perLayer := best.nx * best.ny // number of cells in one full layer
	if perLayer > 0 && best.nLayers > 0 && best.countPlaced > 0 {
		partialCount := best.countPlaced - (best.nLayers-1)*perLayer // items in the last (possibly partial) layer
		if partialCount > 0 && partialCount < perLayer {
			cols := best.nx
			rows := best.ny
			cellW := best.ori.W
			cellH := best.ori.H
			layerL := best.ori.L // last layer thickness along length

			fullRows := partialCount / cols // how many full rows are filled in the last layer
			remCols := partialCount % cols  // how many cells are filled in the next row
			filledRows := fullRows
			if remCols > 0 {
				filledRows++
			}
			freeRows := rows - filledRows // rows above the last filled row

			// Hole to the right of the partially filled row (one-row tall).
			if remCols > 0 {
				holeRightW := (cols - remCols) * cellW
				holeRight := Box{W: holeRightW, H: cellH, L: layerL}
				if holeRight.W > 0 && holeRight.H > 0 && holeRight.L > 0 {
					sub = append(sub, holeRight)
				}
			}

			// Hole above all filled rows (covers full used width).
			if freeRows > 0 {
				holeAboveH := freeRows * cellH
				holeAbove := Box{W: cols * cellW, H: holeAboveH, L: layerL}
				if holeAbove.W > 0 && holeAbove.H > 0 && holeAbove.L > 0 {
					sub = append(sub, holeAbove)
				}
			}
		}
	}

	// Recurse into larger boxes first.
	for i := 0; i < len(sub); i++ { // selection-sort by volume (descending)
		for j := i + 1; j < len(sub); j++ {
			if vol(sub[j]) > vol(sub[i]) { // compare leftover volumes
				sub[i], sub[j] = sub[j], sub[i] // swap to bring larger box first
			}
		}
	}

	for _, b := range sub { // recurse into each leftover sub-box
		ps, rem := packBoxRecursive(b, items)  // recursive packing into leftover region
		placements = append(placements, ps...) // collect placements from recursion
		items = rem                            // carry updated remaining counts forward
	}

	return placements, items // placements made in this box + updated remaining counts
}

type PackSummary struct {
	TotalContainers int            // total number of containers used
	ByContainer     map[string]int // breakdown by container type name
	Remaining       []ItemType     // items left unplaced (if packing becomes impossible)
}

// chooseNextContainer:
// - if remaining items fit into ONE container of some type, choose the SMALLEST such container
// - otherwise choose the container that places the MOST items in one container (tie -> smaller container)
func chooseNextContainer(containers []Container, remaining []ItemType) (Container, []ItemType, int, bool) {
	foundFitAll := false         // whether we found any container that can fit ALL remaining items in one go
	var bestFitAllC Container    // best (smallest) container that fits all remaining items
	var bestFitAllRem []ItemType // remaining items after packing into bestFitAllC (should be all zeros)

	foundPartial := false         // whether we found any container that can place at least something
	var bestPartialC Container    // best container under "partial" rule
	var bestPartialRem []ItemType // remaining items after packing one bestPartialC
	bestPartialPlaced := 0        // how many items bestPartialC places in one container

	for _, c := range containers { // evaluate each container type for the current remaining items
		ps, rem := packBoxRecursive(c.Dim, copyItems(remaining)) // simulate packing ONE container of type c
		placed := sumPlaced(ps)                                  // total items placed into that one container
		if placed == 0 {
			continue
		}

		if allPlaced(rem) { // this container can fit ALL remaining items
			if !foundFitAll || containerVol(c) < containerVol(bestFitAllC) { // choose the smallest such container
				foundFitAll = true
				bestFitAllC = c
				bestFitAllRem = rem
			}
			continue
		}

		if !foundPartial || // first partial candidate
			placed > bestPartialPlaced || // better by placing more items
			(placed == bestPartialPlaced && containerVol(c) < containerVol(bestPartialC)) { // tie-break: smaller container
			foundPartial = true
			bestPartialC = c
			bestPartialRem = rem
			bestPartialPlaced = placed
		}
	}

	if foundFitAll {
		placed := sumPlacedFromCounts(remaining) - sumPlacedFromCounts(bestFitAllRem) // how many were placed overall
		return bestFitAllC, bestFitAllRem, placed, true                               // pick smallest "fits-all" container
	}
	if foundPartial {
		return bestPartialC, bestPartialRem, bestPartialPlaced, true // pick best partial container
	}
	return Container{}, nil, 0, false // nothing can be placed in any container
}

func sumPlacedFromCounts(items []ItemType) int {
	total := 0
	for _, it := range items { // sum remaining counts (used to derive "placed" in fits-all branch)
		total += it.Count
	}
	return total
}

// PackUntilDoneAdaptive keeps packing containers, but if the remaining items can fit into a smaller
// container, it switches to that container for the final container(s).
// Returns total containers used, breakdown by container name, and remaining items (if some cannot be placed).
func PackUntilDoneAdaptive(containers []Container, items []ItemType) PackSummary {
	remaining := copyItems(items) // working copy (we mutate Count during packing)
	s := PackSummary{             // aggregated result
		ByContainer: make(map[string]int), // container usage counters
		Remaining:   remaining,            // will be updated as we pack
	}

	for !allPlaced(remaining) { // keep packing until all counts are zero
		c, rem, placed, ok := chooseNextContainer(containers, remaining) // pick best next container and simulate packing it
		if !ok || placed == 0 {                                          // safety: nothing can be placed, stop to avoid infinite loop
			s.Remaining = remaining
			return s
		}
		s.TotalContainers++     // increment total containers used
		s.ByContainer[c.Name]++ // increment usage for chosen container type
		remaining = rem         // advance to new remaining counts
	}

	s.Remaining = remaining // should be all zeros on success
	return s
}

func main() {
	containers := []Container{ // available container types to choose from
		{Name: "40ft", Dim: Box{W: cmToMm(234.8), H: cmToMm(238.44), L: cmToMm(1203.1)}}, // container option #1
		{Name: "10ft", Dim: Box{W: cmToMm(234.8), H: cmToMm(238.44), L: cmToMm(279.4)}},  // container option #2
	}

	type TestCase struct {
		Name  string     // test case name (printed in output)
		Items []ItemType // list of item types for this test case
	}

	testCases := []TestCase{
		{
			Name: "Case 1: single item type (27 rectangular boxes)",
			Items: []ItemType{
				{Name: "RectangleBox", Count: 27, Dim: Box{
					W: cmToMm(78), // item width in mm
					H: cmToMm(79), // item height in mm
					L: cmToMm(93), // item length in mm
				}},
			},
		},
		{
			Name: "Case 2: mixed items (balls + packages)",
			Items: []ItemType{
				{Name: "ball_R40(cube80)", Count: 50, Dim: Box{W: cmToMm(80), H: cmToMm(80), L: cmToMm(80)}},
				{Name: "pkg_80_100_200", Count: 56, Dim: Box{W: cmToMm(80), H: cmToMm(100), L: cmToMm(200)}},
				{Name: "pkg_60_80_150", Count: 48, Dim: Box{W: cmToMm(60), H: cmToMm(80), L: cmToMm(150)}},
			},
		},
		{
			Name: "Case 3: tiny remainder (2 boxes + 1 small box)",
			Items: []ItemType{
				{Name: "BigBox", Count: 2, Dim: Box{W: cmToMm(200), H: cmToMm(100), L: cmToMm(200)}},
				{Name: "SmallBox", Count: 1, Dim: Box{W: cmToMm(50), H: cmToMm(50), L: cmToMm(50)}},
			},
		},
		{
			Name: "Case 4: oversize item (should not fit in any orientation)",
			Items: []ItemType{
				// One dimension is longer than the longest container length (40ft is ~1203.1cm),
				// so this item should not fit in any container in any orientation.
				{Name: "TooLong", Count: 1, Dim: Box{W: cmToMm(1300), H: cmToMm(10), L: cmToMm(10)}},
			},
		},
		{
			Name: "Case 5: both containers used (40ft then 10ft)",
			Items: []ItemType{
				// For a 80cm cube:
				// - 40ft fits: floor(234.8/80)*floor(238.44/80)*floor(1203.1/80) = 2*2*15 = 60
				// - 10ft fits: floor(234.8/80)*floor(238.44/80)*floor(279.4/80)  = 2*2*3  = 12
				// With 61 cubes: first container will be 40ft (places 60), remaining 1 fits in 10ft.
				{Name: "Cube80", Count: 61, Dim: Box{W: cmToMm(80), H: cmToMm(80), L: cmToMm(80)}},
			},
		},
	}

	for _, tc := range testCases { // run all test cases
		fmt.Printf("\n=== %s ===\n", tc.Name)

		s := PackUntilDoneAdaptive(containers, tc.Items)             // run adaptive packing until all Counts reach 0 (or get stuck)
		fmt.Printf("Total containers used: %d\n", s.TotalContainers) // overall containers count
		fmt.Printf("Breakdown:\n")
		for _, c := range containers { // print breakdown in the same order as input
			if n := s.ByContainer[c.Name]; n > 0 { // skip container types not used
				fmt.Printf("- %s: %d\n", c.Name, n)
			}
		}

		// debug only, when none container can fit any of the boxes
		if !allPlaced(s.Remaining) { // if packing could not complete, print remaining (unplaced) counts
			fmt.Printf("Unplaced (could not fit):\n")
			for _, it := range s.Remaining { // list only types with remaining Count
				if it.Count > 0 { // still unplaced
					fmt.Printf("- %s: %d\n", it.Name, it.Count)
				}
			}
		}
	}
}
