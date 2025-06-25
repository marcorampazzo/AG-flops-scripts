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