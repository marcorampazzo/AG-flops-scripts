import homogeneous_bundles as hb
import exts_on_tilting_bundles as ext

# koszul = hb.HomogeneousDirectSum([])

def euler_form(bundle_1, bundle_2, first_is_torsion=False, second_is_torsion=False):
    if not first_is_torsion and not second_is_torsion:
        expansion = bundle_2 * bundle_1.dual() # TODO: if i reverse the order of the product it gives an error, this shouldn't happen!! Ever!!
        print(expansion)
        return expansion.euler()
    if not first_is_torsion and second_is_torsion:
        if isinstance(bundle_1, hb.HomogeneousIrreducible):
            bundle_1.twist = 0
        else:
            for item in bundle_1:
                item.twist = 0
        if isinstance(bundle_2, hb.HomogeneousIrreducible):
            bundle_2.twist = 0
        else:
            for item in bundle_2:
                item.twist = 0
        return euler_form(bundle_1, bundle_2) # here no torsion, we just restricted
    if first_is_torsion and not second_is_torsion:
        return euler_form(bundle_2, bundle_1 * hb.HomogeneousIrreducible(bundle_1.k, bundle_1.n, [], [4 for i in range(bundle_1.n - bundle_1.k)]), False, True)
    else:
        pass # with the koszul thing




Udual = hb.HomogeneousIrreducible(2, 6, [1], [])

print(euler_form(Udual, Udual, False, True))


# print(Udual.cohomology())

# print((Udual*Udual*Udual).euler())

# FirstWindowG26 = hb.HomogeneousDirectSum([
#     hb.HomogeneousIrreducible(2, 6, [1], [1,1,1,1], 1),
#     hb.HomogeneousIrreducible(2, 6, [], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [1], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [2], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [1,1], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [2,1], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [3,1], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [2,2], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [3,2], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [4,2], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [3,3], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [4,3], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [4,4], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [5,4], [], 1),
#     hb.HomogeneousIrreducible(2, 6, [5,5], [], 1) 
#     ]
# )

# print(ext.is_it_tilting(2, 6, FirstWindowG26, 'U*(-2)', 7))