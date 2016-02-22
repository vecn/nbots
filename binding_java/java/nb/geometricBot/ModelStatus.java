package nb.geometricBot;

public enum ModelStatus {
    OK,
    ZERO_VERTICES,
    ZERO_EDGES,
    REPEATED_VERTICES,
    INCOHERENT_EDGES,
    REPEATED_EDGES, 
    INTERSECTED_EDGES,
    VTX_INTERSECTING_EDGE,
    UNCLOSED_BOUNDARY
}
