using BioSequences.FASTA, DataFrames, CSV, nucpos, Statistics

sib_pos = "/bench/PhD/danpos_results/pooled/sib.Fnor.smooth.positions.xls"
rys_pos = "/bench/PhD/danpos_results/pooled/rys.Fnor.smooth.positions.xls"
sib_refpos = "/bench/PhD/danpos_results/pooled/sib.Fnor.smooth.positions.ref_adjust.xls"
rys_refpos = "/bench/PhD/danpos_results/pooled/rys.Fnor.smooth.positions.ref_adjust.xls"

sib_df=CSV.read(sib_pos)
rys_df=CSV.read(rys_pos)
sib_refdf=CSV.read(sib_refpos)
rys_refdf=CSV.read(rys_refpos)

Threads.@threads for t in 1:4
    t==1 && nucpos.map_positions!(sib_df, rys_df)
    t==2 && nucpos.map_positions!(rys_df, sib_df)
    t==3 && nucpos.map_positions!(sib_refdf, rys_refdf)
    t==4 && nucpos.map_positions!(rys_refdf, sib_refdf)
end

println("Unmapped positions: sib $(length(findall(iszero,sib_df.mapped_pos))), sib_ref $(length(findall(iszero,sib_refdf.mapped_pos)))")
println("Unmapped positions: rys $(length(findall(iszero,rys_df.mapped_pos))), rys_ref $(length(findall(iszero,rys_refdf.mapped_pos)))")
println("Avg overlap: sib $(mean(sib_df.rel_overlap)), sib_ref $(mean(sib_refdf.rel_overlap))")
println("Avg overlap: rys $(mean(rys_df.rel_overlap)), rys_ref $(mean(rys_refdf.rel_overlap))")