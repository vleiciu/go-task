import kotlin.math.roundToInt

// ----------------------------
// Data model (all dimensions in millimeters)
// ----------------------------

data class Box(val w: Int, val h: Int, val l: Int) // mm

data class ItemType(
    val name: String, // item type name (identifier for printing/debugging)
    var count: Int,   // how many items of this type remain to be packed
    val dim: Box      // original item dimensions (in mm)
)

data class Container(
    val name: String, // container type name (identifier for printing/debugging)
    val dim: Box      // container inner dimensions (in mm)
)

// Orientation is an axis-aligned permutation of an item's dimensions.
// It represents how the item is rotated relative to container axes.
data class Orientation(val w: Int, val h: Int, val l: Int)

data class Placement(
    val inBox: Box,          // the sub-box (region) we packed into (size only; no position tracking)
    val item: String,        // item type name that was placed
    val ori: Orientation,    // chosen orientation for this placement
    val count: Int,          // how many items were placed in this placement
    val nx: Int,             // grid count along W in the cross-section
    val ny: Int,             // grid count along H in the cross-section
    val nLayers: Int,        // number of layers along L used by this placement
    val used: Box            // used block size inside inBox (nx*w, ny*h, nLayers*l)
)

data class PackSummary(
    val totalContainers: Int,          // total number of containers used
    val byContainer: Map<String, Int>, // breakdown by container type name
    val remaining: List<ItemType>      // items left unplaced (if packing becomes impossible)
)

// ----------------------------
// Helpers
// ----------------------------

fun cmToMm(x: Double): Int = (x * 10.0).roundToInt() // cm -> mm (rounded)

fun vol(b: Box): Long = b.w.toLong() * b.h.toLong() * b.l.toLong() // volume in mm^3

fun containerVol(c: Container): Long = vol(c.dim) // container volume in mm^3 (used for tie-breaking)

fun ceilDiv(a: Int, b: Int): Int {
    if (b <= 0) return 0
    return (a + b - 1) / b
}

fun minInt(a: Int, b: Int): Int = if (a < b) a else b

fun allPlaced(items: List<ItemType>): Boolean = items.all { it.count <= 0 }

fun sumCounts(items: List<ItemType>): Int = items.sumOf { it.count }

fun sumPlaced(placements: List<Placement>): Int = placements.sumOf { it.count }

fun copyItems(items: List<ItemType>): MutableList<ItemType> =
    items.map { it.copy(count = it.count) }.toMutableList() // deep copy of counts

fun perm3(a: Int, b: Int, c: Int): List<Triple<Int, Int, Int>> = listOf(
    Triple(a, b, c), Triple(a, c, b),
    Triple(b, a, c), Triple(b, c, a),
    Triple(c, a, b), Triple(c, b, a)
)

// ----------------------------
// Minimal JSON parser (no external dependencies)
// ----------------------------

private class JsonParser(private val s: String) {
    private var i: Int = 0 // current cursor position

    fun parse(): Any? {
        skipWs()
        val v = parseValue()
        skipWs()
        if (i != s.length) error("Trailing data at position $i")
        return v
    }

    private fun parseValue(): Any? {
        skipWs()
        if (i >= s.length) error("Unexpected end of input")
        return when (val ch = s[i]) {
            '{' -> parseObject()
            '[' -> parseArray()
            '"' -> parseString()
            't' -> { expect("true"); true }
            'f' -> { expect("false"); false }
            'n' -> { expect("null"); null }
            '-', in '0'..'9' -> parseNumber()
            else -> error("Unexpected character '$ch' at position $i")
        }
    }

    private fun parseObject(): Map<String, Any?> {
        expectChar('{')
        skipWs()
        val out = linkedMapOf<String, Any?>()
        if (peekChar('}')) {
            expectChar('}')
            return out
        }
        while (true) {
            skipWs()
            val key = parseString()
            skipWs()
            expectChar(':')
            val value = parseValue()
            out[key] = value
            skipWs()
            if (peekChar('}')) {
                expectChar('}')
                return out
            }
            expectChar(',')
        }
    }

    private fun parseArray(): List<Any?> {
        expectChar('[')
        skipWs()
        val out = arrayListOf<Any?>()
        if (peekChar(']')) {
            expectChar(']')
            return out
        }
        while (true) {
            val value = parseValue()
            out.add(value)
            skipWs()
            if (peekChar(']')) {
                expectChar(']')
                return out
            }
            expectChar(',')
        }
    }

    private fun parseString(): String {
        expectChar('"')
        val sb = StringBuilder()
        while (i < s.length) {
            val ch = s[i++]
            when (ch) {
                '"' -> return sb.toString()
                '\\' -> {
                    if (i >= s.length) error("Unterminated escape at position $i")
                    val esc = s[i++]
                    when (esc) {
                        '"', '\\', '/' -> sb.append(esc)
                        'b' -> sb.append('\b')
                        'f' -> sb.append('\u000C')
                        'n' -> sb.append('\n')
                        'r' -> sb.append('\r')
                        't' -> sb.append('\t')
                        'u' -> {
                            if (i + 4 > s.length) error("Invalid unicode escape at position $i")
                            val hex = s.substring(i, i + 4)
                            i += 4
                            sb.append(hex.toInt(16).toChar())
                        }
                        else -> error("Invalid escape '\\$esc' at position ${i - 1}")
                    }
                }
                else -> sb.append(ch)
            }
        }
        error("Unterminated string")
    }

    private fun parseNumber(): Double {
        val start = i
        if (s[i] == '-') i++
        while (i < s.length && s[i].isDigit()) i++
        if (i < s.length && s[i] == '.') {
            i++
            while (i < s.length && s[i].isDigit()) i++
        }
        if (i < s.length && (s[i] == 'e' || s[i] == 'E')) {
            i++
            if (i < s.length && (s[i] == '+' || s[i] == '-')) i++
            while (i < s.length && s[i].isDigit()) i++
        }
        return s.substring(start, i).toDouble()
    }

    private fun skipWs() {
        while (i < s.length && s[i].isWhitespace()) i++
    }

    private fun expectChar(c: Char) {
        if (i >= s.length || s[i] != c) error("Expected '$c' at position $i")
        i++
    }

    private fun peekChar(c: Char): Boolean = i < s.length && s[i] == c

    private fun expect(word: String) {
        if (!s.regionMatches(i, word, 0, word.length)) error("Expected '$word' at position $i")
        i += word.length
    }
}

private fun Any?.asObj(): Map<String, Any?> = this as? Map<String, Any?> ?: error("Expected object")
private fun Any?.asArr(): List<Any?> = this as? List<Any?> ?: error("Expected array")
private fun Any?.asStr(): String = this as? String ?: error("Expected string")
private fun Any?.asNum(): Double = when (this) {
    is Double -> this
    is Int -> this.toDouble()
    is Long -> this.toDouble()
    else -> error("Expected number")
}

// ----------------------------
// AllOrientations
// ----------------------------

// AllOrientations returns all unique orientations of the item dimensions,
// aligned to container axes (W, H, L).
fun allOrientations(it: ItemType): List<Orientation> {
    val seen = HashSet<Triple<Int, Int, Int>>() // set of (W,H,L) to remove duplicates (e.g., cubes)
    val out = ArrayList<Orientation>(6)         // up to 6 unique permutations

    for (p in perm3(it.dim.w, it.dim.h, it.dim.l)) { // generate all 3D permutations
        val o = Orientation(p.first, p.second, p.third) // treat permuted dims as axis-aligned orientation
        val key = Triple(o.w, o.h, o.l)                 // de-duplication key
        if (seen.add(key)) {
            out.add(o) // collect unique orientation
        }
    }
    return out // all unique axis-aligned orientations
}

// ----------------------------
// Recursive packing into a box (greedy + guillotine leftovers)
// ----------------------------

data class BestCandidate(
    val itemIdx: Int,         // index of chosen item type in items list
    val ori: Orientation,     // chosen orientation for this candidate
    val countPlaced: Int,     // how many items we can place (bounded by Count and capacity)
    val nx: Int,              // grid fit along W
    val ny: Int,              // grid fit along H
    val nLayers: Int,         // layers used along L
    val used: Box,            // block size used by this candidate (in mm)
    val unusedVol: Long       // heuristic: (box volume - placed*item volume)
)

fun betterCandidate(a: BestCandidate, b: BestCandidate): Boolean {
    if (a.countPlaced != b.countPlaced) return a.countPlaced > b.countPlaced // maximize placed items
    if (a.unusedVol != b.unusedVol) return a.unusedVol < b.unusedVol         // minimize leftover volume
    val ua = vol(a.used)
    val ub = vol(b.used)
    if (ua != ub) return ua > ub                                              // prefer more compact usage
    return a.nLayers > b.nLayers                                              // prefer more layers as a tie-breaker
}

// packBoxRecursive greedily packs items into 'box' using axis-aligned grid placements,
// then recursively tries to pack remaining items into leftover sub-boxes.
//
// Leftover partition is a simple non-overlapping guillotine split:
// - right: leftover along W for the used length region
// - top: leftover along H for the used length region (only over usedW)
// - back: leftover along L for the full cross-section
//
// NOTE: This function MUTATES 'items' counts in-place and returns placements done in this box.
fun packBoxRecursive(box: Box, items: MutableList<ItemType>): List<Placement> {
    var found = false                   // whether any candidate placement was found for this box
    var best: BestCandidate? = null     // best candidate found so far

    for ((i, it) in items.withIndex()) { // iterate over item types
        if (it.count <= 0) continue      // skip items with nothing left to pack

        for (o in allOrientations(it)) { // try all axis-aligned orientations
            if (o.w <= 0 || o.h <= 0 || o.l <= 0) continue // guard against invalid dimensions
            if (o.w > box.w || o.h > box.h || o.l > box.l) continue // must fit within current box

            val nx = box.w / o.w // how many items fit along width
            val ny = box.h / o.h // how many items fit along height
            val nz = box.l / o.l // how many items fit along length (max layers)
            if (nx <= 0 || ny <= 0 || nz <= 0) continue

            val capacity = nx * ny * nz // maximum items for this orientation in this box
            val placed = minInt(it.count, capacity) // actual items placed limited by remaining Count
            if (placed <= 0) continue

            val perLayer = nx * ny // items per single layer along length
            val nLayers = minInt(nz, ceilDiv(placed, perLayer)) // layers actually needed to place 'placed' items
            if (nLayers <= 0) continue

            val used = Box(w = nx * o.w, h = ny * o.h, l = nLayers * o.l) // bounding block we occupy
            val itemVol = vol(Box(o.w, o.h, o.l))
            var unusedVol = vol(box) - placed.toLong() * itemVol // heuristic leftover volume
            if (unusedVol < 0) unusedVol = 0

            val cand = BestCandidate(
                itemIdx = i,
                ori = o,
                countPlaced = placed,
                nx = nx,
                ny = ny,
                nLayers = nLayers,
                used = used,
                unusedVol = unusedVol
            )

            if (!found || betterCandidate(cand, best!!)) {
                found = true
                best = cand
            }
        }
    }

    if (!found || best == null) return emptyList() // nothing fits into this box

    // Apply placement: decrement count.
    val chosen = best!!
    items[chosen.itemIdx].count -= chosen.countPlaced // update remaining Count for the chosen item type
    if (items[chosen.itemIdx].count < 0) items[chosen.itemIdx].count = 0

    val placements = ArrayList<Placement>()
    placements.add(
        Placement(
            inBox = box,                   // box we packed into
            item = items[chosen.itemIdx].name, // item type name
            ori = chosen.ori,              // orientation chosen
            count = chosen.countPlaced,    // how many were placed
            nx = chosen.nx,                // grid along W
            ny = chosen.ny,                // grid along H
            nLayers = chosen.nLayers,      // layers along L
            used = chosen.used             // occupied block dimensions
        )
    )

    // Build leftover sub-boxes (non-overlapping).
    val usedW = chosen.used.w // occupied extent along W for splitting
    val usedH = chosen.used.h // occupied extent along H for splitting
    val usedL = chosen.used.l // occupied extent along L for splitting

    val right = Box(w = box.w - usedW, h = box.h, l = usedL) // leftover region to the "right" (along W)
    val top = Box(w = usedW, h = box.h - usedH, l = usedL)   // leftover region on "top" (along H)
    val back = Box(w = box.w, h = box.h, l = box.l - usedL)  // leftover region "behind" (along L)

    val sub = ArrayList<Box>(3) // list of leftover sub-boxes we will recurse into
    if (right.w > 0 && right.h > 0 && right.l > 0) sub.add(right)
    if (top.w > 0 && top.h > 0 && top.l > 0) sub.add(top)
    if (back.w > 0 && back.h > 0 && back.l > 0) sub.add(back)

    // If the last layer is only partially filled, expose the internal empty cells of that last layer
    // as additional leftover sub-boxes (so we can "finish" the layer with other items).
    //
    // IMPORTANT: We do not track coordinates, only sizes. These sub-boxes are still non-overlapping
    // with right/top/back because they are inside the used cross-section (usedW×usedH) and only
    // within the last layer thickness (chosen.ori.l).
    val perLayer = chosen.nx * chosen.ny // number of cells in one full layer
    if (perLayer > 0 && chosen.nLayers > 0 && chosen.countPlaced > 0) {
        val partialCount = chosen.countPlaced - (chosen.nLayers - 1) * perLayer // items in the last (possibly partial) layer
        if (partialCount > 0 && partialCount < perLayer) {
            val cols = chosen.nx
            val rows = chosen.ny
            val cellW = chosen.ori.w
            val cellH = chosen.ori.h
            val layerL = chosen.ori.l // last layer thickness along length

            val fullRows = partialCount / cols // how many full rows are filled in the last layer
            val remCols = partialCount % cols  // how many cells are filled in the next row
            var filledRows = fullRows
            if (remCols > 0) filledRows++
            val freeRows = rows - filledRows // rows above the last filled row

            // Hole to the right of the partially filled row (one-row tall).
            if (remCols > 0) {
                val holeRightW = (cols - remCols) * cellW
                val holeRight = Box(w = holeRightW, h = cellH, l = layerL)
                if (holeRight.w > 0 && holeRight.h > 0 && holeRight.l > 0) sub.add(holeRight)
            }

            // Hole above all filled rows (covers full used width).
            if (freeRows > 0) {
                val holeAboveH = freeRows * cellH
                val holeAbove = Box(w = cols * cellW, h = holeAboveH, l = layerL)
                if (holeAbove.w > 0 && holeAbove.h > 0 && holeAbove.l > 0) sub.add(holeAbove)
            }
        }
    }

    // Recurse into larger boxes first (by volume, descending).
    sub.sortByDescending { vol(it) }

    for (b in sub) { // recurse into each leftover sub-box
        val ps = packBoxRecursive(b, items) // recursive packing into leftover region
        placements.addAll(ps)               // collect placements from recursion
    }

    return placements
}

// ----------------------------
// Adaptive container choice + packing until all counts reach 0
// ----------------------------

// chooseNextContainer:
// - if remaining items fit into ONE container of some type, choose the SMALLEST such container
// - otherwise choose the container that places the MOST items in one container (tie -> smaller container)
fun chooseNextContainer(
    containers: List<Container>,
    remaining: List<ItemType>
): Triple<Container, MutableList<ItemType>, Int>? {
    var foundFitAll = false
    var bestFitAllC: Container? = null
    var bestFitAllRem: MutableList<ItemType>? = null

    var foundPartial = false
    var bestPartialC: Container? = null
    var bestPartialRem: MutableList<ItemType>? = null
    var bestPartialPlaced = 0

    for (c in containers) { // evaluate each container type for the current remaining items
        val tmp = copyItems(remaining) // simulate on a fresh copy so we don't mutate the caller
        val ps = packBoxRecursive(c.dim, tmp) // simulate packing ONE container of type c
        val placed = sumPlaced(ps) // total items placed into that one container
        if (placed == 0) continue

        if (allPlaced(tmp)) { // this container can fit ALL remaining items
            if (!foundFitAll || containerVol(c) < containerVol(bestFitAllC!!)) {
                foundFitAll = true
                bestFitAllC = c
                bestFitAllRem = tmp
            }
            continue
        }

        if (!foundPartial ||
            placed > bestPartialPlaced ||
            (placed == bestPartialPlaced && containerVol(c) < containerVol(bestPartialC!!))
        ) {
            foundPartial = true
            bestPartialC = c
            bestPartialRem = tmp
            bestPartialPlaced = placed
        }
    }

    return when {
        // Pick smallest "fits-all" container; compute how many items got placed in this one container.
        foundFitAll -> Triple(
            bestFitAllC!!,
            bestFitAllRem!!,
            sumCounts(remaining) - sumCounts(bestFitAllRem!!)
        )
        // Pick best partial container; reuse the already computed bestPartialPlaced.
        foundPartial -> Triple(bestPartialC!!, bestPartialRem!!, bestPartialPlaced)
        else -> null // nothing can be placed in any container
    }
}

// PackUntilDoneAdaptive keeps packing containers, but if the remaining items can fit into a smaller
// container, it switches to that container for the final container(s).
fun packUntilDoneAdaptive(containers: List<Container>, items: List<ItemType>): PackSummary {
    var remaining = copyItems(items) // working copy (we mutate Count during packing)
    val byContainer = linkedMapOf<String, Int>() // container usage counters (stable order)
    var totalContainers = 0

    while (!allPlaced(remaining)) { // keep packing until all counts are zero
        val choice = chooseNextContainer(containers, remaining) ?: break // pick best next container and simulate packing it
        val c = choice.first   // chosen container type for this step
        val rem = choice.second // remaining counts after packing ONE container of type c
        val placed = choice.third // how many items were placed into that one container

        // Safety: if nothing changes (shouldn't happen), stop to avoid infinite loop.
        if (placed <= 0) break

        totalContainers++
        byContainer[c.name] = (byContainer[c.name] ?: 0) + 1
        remaining = rem
    }

    return PackSummary(
        totalContainers = totalContainers,
        byContainer = byContainer,
        remaining = remaining
    )
}

// ----------------------------
// Demo input / output (test cases loaded from JSON)
// ----------------------------

data class TestCase(
    val name: String,         // test case name (printed in output)
    val items: List<ItemType> // list of item types for this test case
)

val inputPath = if (args.isNotEmpty()) args[0] else "testcases.json" // JSON file path (default: repo file)
val rawJson = java.io.File(inputPath).readText() // read JSON content from disk
val root = JsonParser(rawJson).parse().asObj()   // parse JSON into a generic object tree

val containers = root["containers"].asArr().map { cAny ->
    val c = cAny.asObj()
    val name = c["name"].asStr()
    val dim = c["dimCm"].asObj()
    Container(
        name = name,
        dim = Box(
            w = cmToMm(dim["w"].asNum()),
            h = cmToMm(dim["h"].asNum()),
            l = cmToMm(dim["l"].asNum())
        )
    )
}

val testCases = root["cases"].asArr().map { tcAny ->
    val tc = tcAny.asObj()
    val name = tc["name"].asStr()
    val items = tc["items"].asArr().map { itAny ->
        val it = itAny.asObj()
        val dim = it["dimCm"].asObj()
        ItemType(
            name = it["name"].asStr(),
            count = it["count"].asNum().toInt(),
            dim = Box(
                w = cmToMm(dim["w"].asNum()),
                h = cmToMm(dim["h"].asNum()),
                l = cmToMm(dim["l"].asNum())
            )
        )
    }
    TestCase(name = name, items = items)
}

for (tc in testCases) {
    println("\n=== ${tc.name} ===")

    val summary = packUntilDoneAdaptive(containers, tc.items) // run adaptive packing until all Counts reach 0 (or get stuck)

    println("Total containers used: ${summary.totalContainers}") // overall containers count
    println("Breakdown:")
    for (c in containers) { // print breakdown in the same order as input
        val n = summary.byContainer[c.name] ?: 0
        if (n > 0) println("- ${c.name}: $n")
    }

    // Debug only, when introduced boxes cannot fit in any container
    if (!allPlaced(summary.remaining)) { // if packing could not complete, print remaining (unplaced) counts
        println("Unplaced (could not fit):")
        for (it in summary.remaining) { // list only types with remaining Count
            if (it.count > 0) println("- ${it.name}: ${it.count}")
        }
    }
}


