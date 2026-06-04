import homogeneous_bundles as hb
import exts_on_tilting_bundles as ext

koszul_complex = hb.HomogeneousDirectSum([
    hb.HomogeneousIrreducible(2, 6, [], [], 1, 0, -1),
    hb.HomogeneousIrreducible(2, 6, [], [1], 1, -1, -2),
    hb.HomogeneousIrreducible(2, 6, [], [1, 1], 1, -2, -3),
    hb.HomogeneousIrreducible(2, 6, [], [1, 1, 1], 1, -3, -4),
    hb.HomogeneousIrreducible(2, 6, [], [1, 1, 1, 1], 1, -4, -5),
])

def euler_form(bundle_1, bundle_2, first_is_torsion=False, second_is_torsion=False):
    if not first_is_torsion and not second_is_torsion:
        expansion = bundle_2 * bundle_1.dual() # TODO: if i reverse the order of the product it gives an error, this shouldn't happen!! Ever!!
        # print(expansion)
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
        bundle_1 = bundle_1 * koszul_complex
        # bundle_2 = bundle_2 * koszul_complex
        # print(bundle_1)
        return euler_form(bundle_1, bundle_2, False, True)



def euler_matrix(list_of_bundles_1, list_of_bundles_2):
    matrix = [
        [euler_form(b1[0], b2[0], b1[1], b2[1]) for b2 in list_of_bundles_2]
        for b1 in list_of_bundles_1
    ]
    # Pretty print
    width = max(len(str(x)) for row in matrix for x in row)
    for row in matrix:
        print(" ".join(f"{x:>{width}}" for x in row))
    return matrix    



# in the following, ther boolean variable indicates whether the bundle is a pushforward from
# the zero section (True) or not (False)

collection = [
    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 0, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 0, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 1, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 1, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 2, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 2, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 3, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 3, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 4, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 4, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 5, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 5, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 6, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 6, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 7, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 7, 0), False],

    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 8, 0), False],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 8, 0), False],
    
    [hb.HomogeneousIrreducible(2, 6, [], [], 1, 0, 0), True],
    [hb.HomogeneousIrreducible(2, 6, [1], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [2], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [1,1], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [2,1], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [3,1], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [2,2], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [3,2], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [4,2], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [3,3], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [4,3], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [4,4], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [5,4], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [5,5], [], 1, 0, 0),True],
    [hb.HomogeneousIrreducible(2, 6, [6,5], [], 1, 0, 0),True],

    
]

euler_matrix(collection, collection)


grassmannian_collection = [
    hb.HomogeneousIrreducible(2, 6, [], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [1], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [2], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [1,1], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [2,1], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [3,1], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [2,2], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [3,2], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [4,2], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [3,3], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [4,3], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [4,4], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [5,4], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [5,5], [], 1, 0, 0),
    hb.HomogeneousIrreducible(2, 6, [6,5], [], 1, 0, 0),
]

# euler_matrix(grassmannian_collection, grassmannian_collection)

# print(euler_form(O, O))
# print(euler_form(O, Udual))
# print(euler_form(Udual, O))
# print(euler_form(Udual, Udual))


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