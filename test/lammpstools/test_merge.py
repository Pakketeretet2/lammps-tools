import lammpstools, dumpreader, block_data

d = dumpreader.dumpreader_cpp( "melt.dump" )
b = d.getblock()
print("Block at t = ", b.meta.t, " has ", b.meta.N, " particles and atom style",
        b.meta.atom_style)

for b2 in d:
    pass

print("Block at t = ", b2.meta.t, " has ", b2.meta.N, " particles and atom style",
        b2.meta.atom_style)
b2.merge(b, error_on_double_id = False)
print("Merging... After merge,")
print("Block at t = ", b2.meta.t, " has ", b2.meta.N, " particles and atom style",
        b2.meta.atom_style)
block_data.block_to_data( b2, "test.data", overwrite = True )

