import homogeneous_bundles as hb
import exts_on_tilting_bundles as ext

# def euler_char(bundle_1, bundle_2):
#     if isinstance(bundle_2, torsion):
#         if isinstance(bundle_1, torsion):
#             pass # Euler char for two torsion bundles
#         else:
#             pass # if only the second is torsion
#     else:
#         pass # dualize and apply koszul to the first one


Udual = hb.HomogeneousIrreducible(2, 6, [1], [])

print(Udual.coh())

print((Udual).euler())

FirstWindowG26 = hb.HomogeneousDirectSum([
    hb.HomogeneousIrreducible(2, 6, [1], [1,1,1,1], 1),
    hb.HomogeneousIrreducible(2, 6, [], [], 1),
    hb.HomogeneousIrreducible(2, 6, [1], [], 1),
    hb.HomogeneousIrreducible(2, 6, [2], [], 1),
    hb.HomogeneousIrreducible(2, 6, [1,1], [], 1),
    hb.HomogeneousIrreducible(2, 6, [2,1], [], 1),
    hb.HomogeneousIrreducible(2, 6, [3,1], [], 1),
    hb.HomogeneousIrreducible(2, 6, [2,2], [], 1),
    hb.HomogeneousIrreducible(2, 6, [3,2], [], 1),
    hb.HomogeneousIrreducible(2, 6, [4,2], [], 1),
    hb.HomogeneousIrreducible(2, 6, [3,3], [], 1),
    hb.HomogeneousIrreducible(2, 6, [4,3], [], 1),
    hb.HomogeneousIrreducible(2, 6, [4,4], [], 1),
    hb.HomogeneousIrreducible(2, 6, [5,4], [], 1),
    hb.HomogeneousIrreducible(2, 6, [5,5], [], 1) 
    ]
)

print(ext.is_it_tilting(2, 6, FirstWindowG26, 'U*(-2)', 7))