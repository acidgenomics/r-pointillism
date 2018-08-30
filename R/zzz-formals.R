# # Set the formals.
# f <- formals(.plotDot)
# fshared <- formals(.plotGene)[
#     intersect(
#         formalArgs(.plotDot),
#         formalArgs(.plotGene)
#     )]
# f[names(fshared)] <- fshared
# formals(.plotDot) <- as.pairlist(f)



# # Set the formals.
# f <- formals(.plotViolin)
# fshared <- formals(.plotGene)[
#     intersect(
#         formalArgs(.plotViolin),
#         formalArgs(.plotGene)
#     )]
# f[names(fshared)] <- fshared
# formals(.plotViolin) <- as.pairlist(f)



# # Set the formals.
# f <- formals(.plotFeature)
# fshared <- formals(.plotMarker)[
#     intersect(
#         formalArgs(.plotFeature),
#         formalArgs(.plotMarker)
#     )]
# f[names(fshared)] <- fshared
# formals(.plotFeature) <- as.pairlist(f)
