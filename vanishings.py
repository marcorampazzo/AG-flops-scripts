import homogeneous_bundles as hb
import exts_on_tilting_bundles as ext


for k in range(5):
    line_bundle = hb.HomogeneousIrreducible(2, 5, [0], [k,k,k])
    print(f"O({-k}) has no higher cohomology: {ext.higher_cohomlogy_on_total_space(2, 5, line_bundle, 'Q(-2)')}")

for k in range(3):
    vector_bundle = hb.HomogeneousIrreducible(2, 5, [1], [k,k,k])
    print(f"U^*({-k}) has no higher cohomology: {ext.higher_cohomlogy_on_total_space(2, 5, vector_bundle, 'Q(-2)')}")


for k in range(1, 4):
    vector_bundle = hb.HomogeneousIrreducible(2, 5, [2+k,2+k,k], [2,2])
    print(f"Sym^2 U_3({k}) has no higher cohomology: {ext.higher_cohomlogy_on_total_space(2, 5, vector_bundle, 'Q(-2)')}")

print(f"Cohomology of Sym^2 U_3(-2) on G(3, 5): {hb.HomogeneousIrreducible(3, 5, [2,2], [4,4]).cohomology()}")



for a in range(1, 3):
    twisted_sym = hb.HomogeneousIrreducible(3, 5, [a-1,a-1], [3,3])
    universal = hb.HomogeneousIrreducible(3, 5, [1], [])
    print(f"Cohomology of U_3^* \otimes Sym^({a-1}) U_3({a-4}) on G(3, 5): {(universal*twisted_sym).cohomology()}")


for a in range(3):
    twisted_sym = hb.HomogeneousIrreducible(2, 5, [2+a,a], [])
    universal = hb.HomogeneousIrreducible(2, 5, [1], [])
    print(f"vanishing of higher Ext(U_2^*, Sym^2 U_2^*({a})) on G(2, 5): {ext.higher_ext_on_total_space(2, 5, universal, twisted_sym, 'Q(-2)')}")    


for a in range(1,3):
    sym = hb.HomogeneousIrreducible(3, 5, [2], [])
    twisted_sym = hb.HomogeneousIrreducible(3, 5, [a-1, a-1], [3,3])
    print(f"cohomology of Sym^2 U_2^* \otimes Sym^{a-1} U_2^*({a-4})) on G(3, 5) = {(sym*twisted_sym).cohomology()}")    

for a in range(2,4):
    universal = hb.HomogeneousIrreducible(3, 5, [1], [0])
    twisted_sym = hb.HomogeneousIrreducible(3, 5, [2+a,2+a, a], [2, 2])
    print(f"vanishing of higher Ext(U_3^*, Sym^2 U_3({a})) on G(3, 5): {ext.higher_ext_on_total_space(3, 5, universal, twisted_sym, 'U*(-2)')}")   


sym = hb.HomogeneousIrreducible(2, 5, [2], [])
dual_sym = hb.HomogeneousIrreducible(2, 5, [2], [2,2,2])
print(f"Sym^2 U_2\otimes \Sym^2 U_2^* has no higher cohomology: {ext.higher_cohomlogy_on_total_space(3, 5, vector_bundle, 'U*(-2)')}")



for a in range(3):
    twisted_sym = hb.HomogeneousIrreducible(3, 5, [2+a,2+a, a], [2, 2])
    universal = hb.HomogeneousIrreducible(3, 5, [1,1], [1,1])
    print(f"vanishing of higher Ext(U_3^*, Sym^2 U_3^*({a})) on G(3, 5): {ext.higher_ext_on_total_space(3, 5, twisted_sym, universal, 'U*(-2)')}")   


sym = hb.HomogeneousIrreducible(3, 5, [2], [])
universal = hb.HomogeneousIrreducible(3, 5, [1], [])
line_bundle = hb.HomogeneousIrreducible(3, 5, [], [3,3])
print(f"Cohomology of U_3^* \otimes Sym^2 U_3^*\otimes O(-3) on G(3, 5): {(sym*universal*line_bundle).cohomology()}")
