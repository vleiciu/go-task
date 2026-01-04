# go-task

## How to run

### Run Go (`main.go`)

**Requirements**
- Go (see `go.mod` for the configured Go version)

**Command**
```bash
cd /Users/vleiciu/Desktop/gotask
go run main.go
```

### Run Kotlin script (`pack.kts`)

**Requirements**
- Kotlin compiler (`kotlinc`)

**Command**
```bash
cd /Users/vleiciu/Desktop/gotask
kotlinc -script pack.kts
```

## What was implemented

This repo contains a simplified packing simulation that estimates how many containers are needed to pack a list of items.

### Key pieces

- **All item orientations**
  - `AllOrientations(it ItemType) []Orientation` (Go) / `allOrientations(it: ItemType): List<Orientation>` (Kotlin)
  - Generates all unique axis-aligned permutations of the item dimensions (up to 6).

- **Greedy recursive packing into leftover space**
  - `packBoxRecursive(box, items)` tries to place items into a box using a simple grid fit:
    - \(nx = \lfloor W_{box}/W_{item}\rfloor\), \(ny = \lfloor H_{box}/H_{item}\rfloor\), \(nz = \lfloor L_{box}/L_{item}\rfloor\)
    - `placed = min(count, nx*ny*nz)`
    - Uses only as many layers along length as needed for `placed`, then splits the remaining space into 3 non-overlapping sub-boxes:
      - `right` (leftover width), `top` (leftover height), `back` (leftover length)
    - Recursively repeats the same process for those sub-boxes.

- **Adaptive container selection + total container counting**
  - `PackUntilDoneAdaptive(containers, items)` repeats “pack ONE container” until all `Count == 0`.
  - At each step it chooses the next container type with `chooseNextContainer(...)`:
    - If the entire remaining set fits into a single container of some type, it picks the **smallest** such container.
    - Otherwise it picks the container that places the **most items** in one container (tie-break: smaller container volume).

### Notes / limitations

- This is a heuristic (not an optimal 3D bin-packing solver).
- Placement is axis-aligned; there is no tracking of exact 3D coordinates—only sizes and counts.