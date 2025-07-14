import homogeneous_bundles as hb

def ext_on_grassmannian(E, F):
    """
    outputs a list of dictionaries of this form:
    [{'representaiton': [0, 0, 0, 0, 0], 'dimension': 1, 'degree': 1}, ....]
    """
    return (E.dual()*F).cohomology()


def ext_on_total_space(k, n, E, F, bundle, cutoff=8):
    """
    compute the Ext of two HomogeneousIrreducible's (or direct sums of them).
    We must input k, n as in G(k, n), the two bundles, the type of total space
    and a cutoff integer which is 8 by default.

    The type of total space is quite restrictive so far, one could make this
    more versatile by replacing the append's with symmetric plethysms, but we don't
    need it so far.

    The cutoff consists in the highest symmetric power we consider in the expansion
    of the pushforward of the structure sheaf wrt the bundle map. After a given power,
    everything becomes globally generated, and for our applications 7 or 8 is more than enough.
    """
    pushforward_structure_sheaf = []
    if bundle == 'U*(-2)':
        for i in range(cutoff + 1):
            pushforward_structure_sheaf.append(hb.HomogeneousIrreducible(k, n, [2*i]*(k-1)+[i], [0]*(n-k), 1))
        pushforward_structure_sheaf = hb.HomogeneousDirectSum(pushforward_structure_sheaf)
        # return pushforward_structure_sheaf
    elif bundle == 'Q(-2)':
        for i in range(cutoff + 1):
            pushforward_structure_sheaf.append(hb.HomogeneousIrreducible(k, n, [2*i]*(k), [i]+[0]*(n-k-1), 1))
        pushforward_structure_sheaf = hb.HomogeneousDirectSum(pushforward_structure_sheaf)
        # return pushforward_structure_sheaf
    else:
        return "ERROR: this bundle is not supported (yet)"
    return (E.dual() * F * pushforward_structure_sheaf).cohomology()

def higher_ext_on_total_space(k, n, E, F, bundle, cutoff=8):
    """
    returns True if there are NO higher Exts.
    """
    check = True
    cohomology = ext_on_total_space(k, n, E, F, bundle, cutoff)
    for item in cohomology:
        if item != "no cohomology":
            if item["degree"] != 0:
                check = False
    return check 

def higher_cohomlogy_on_total_space(k, n, E, bundle, cutoff=8):
    structure_sheaf = hb.HomogeneousIrreducible(k, n, [], [])
    return higher_ext_on_total_space(k, n, structure_sheaf, E, bundle, cutoff) 

def is_it_tilting(k, n, E, bundle, cutoff):
    """
    this function checks if a bundle is partially tilting by
    computing the self Ext on the total space. It gives True if it is.
    """
    check = True
    cohomology = ext_on_total_space(k, n, E, E, bundle, cutoff)
    for item in cohomology:
        if item != "no cohomology":
            if item["degree"] != 0:
                check = False
    return check


# here are our partially tilting bundles

FirstWindowG35 = hb.HomogeneousDirectSum([
    hb.HomogeneousIrreducible(3, 5, [], [], 1),
    hb.HomogeneousIrreducible(3, 5, [1,1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [], [1,1], 1),
    hb.HomogeneousIrreducible(3, 5, [1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [1], [1,1], 1),
    hb.HomogeneousIrreducible(3, 5, [2], [1,1], 1),
    hb.HomogeneousIrreducible(3, 5, [1,1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,2,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,2,2], [], 1),
    ]
)

SecondWindowG35 = hb.HomogeneousDirectSum([
    hb.HomogeneousIrreducible(3, 5, [], [1,1], 1),
    hb.HomogeneousIrreducible(3, 5, [], [], 1),
    hb.HomogeneousIrreducible(3, 5, [1,1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,2,2], [], 1),
    hb.HomogeneousIrreducible(3, 5, [3,3,3], [], 1),
    hb.HomogeneousIrreducible(3, 5, [1], [1,1,1], 1),
    hb.HomogeneousIrreducible(3, 5, [1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [1,1], [], 1),
    hb.HomogeneousIrreducible(3, 5, [2,2,1], [1,1], 1),
    ]
)

# and here we check that these bundles are indeed partially tilting
if __name__ == '__main__':
    print(f"Is the first generator tilting? {is_it_tilting(3, 5, FirstWindowG35, 'U*(-2)', 7)}.")
    print(f"Is the second generator tilting? {is_it_tilting(3, 5, SecondWindowG35, 'U*(-2)', 7)}.")


