@setup_workload begin

    mp2d_ξ   = rand(2, 20)
    mp3d_ξ   = rand(3, 20)
    ix       = 1
    dev_name = :cpu
    df1, df2, df3, df4, df5, df6, df7, df8, df9, T2_1 = rand(10)
    
    @compile_workload begin
        quiet() do
            # utils.jl
            format_seconds(12345)
            
            # type.jl - test both precision types
            for precision in [:single, :double]
                ϵ = set_precision(precision)
                
                # test all dimensions and basis functions
                for (dim, bases) in [(Dim2D(), [Linear(), uGIMP(), Bspline(), Cspline()]),
                                     (Dim3D(), [Linear(), uGIMP(), Bspline(), Cspline()])]
                    for basis in bases
                        # HDF5 configurations
                        h5_no = noHDF5()
                        h5_yes = hasHDF5(interval=10, varnames=(:F, :ξ, :σij, :vs))
                        
                        for h5 in [h5_no, h5_yes]
                            conf = generate_config(ϵ, dim, basis, dev_name, h5, t_tol=1.0, Δt=0.1, 
                                prjpath=tempdir(), prjname="MPS_Precompile_$(precision)_$(typeof(dim))_$(typeof(basis))")
                        end
                    end
                end
                
                # generate field structures with different data types
                if precision === :double
                    mp2d = generate_fields(ϵ, F=rand(4, 20), ξ=mp2d_ξ, σij=rand(4, 20), σm=rand(20),
                        vs=rand(2, 20), ρs=rand(20), Ω=rand(20), nid=ones(Int64, 20))
                    mp3d = generate_fields(ϵ, F=rand(9, 20), ξ=mp3d_ξ, σij=rand(6, 20), σm=rand(20),
                        vs=rand(3, 20), ρs=rand(20), Ω=rand(20), nid=ones(Int64, 20))
                else
                    mp2d = generate_fields(ϵ, F=rand(Float32, 4, 20), ξ=Float32.(mp2d_ξ), σij=rand(Float32, 4, 20), 
                        σm=rand(Float32, 20), vs=rand(Float32, 2, 20), ρs=rand(Float32, 20), 
                        Ω=rand(Float32, 20), nid=ones(Int32, 20))
                    mp3d = generate_fields(ϵ, F=rand(Float32, 9, 20), ξ=Float32.(mp3d_ξ), σij=rand(Float32, 6, 20),
                        σm=rand(Float32, 20), vs=rand(Float32, 3, 20), ρs=rand(Float32, 20),
                        Ω=rand(Float32, 20), nid=ones(Int32, 20))
                end
                
                grid2d = generate_fields(ϵ, 1:10, 1:10, vs=rand(precision === :double ? Float64 : Float32, 2, 100),
                    ms=rand(precision === :double ? Float64 : Float32, 100))
                grid3d = generate_fields(ϵ, 1:5, 1:5, 1:5, vs=rand(precision === :double ? Float64 : Float32, 3, 125),
                    ms=rand(precision === :double ? Float64 : Float32, 125))
                attr = generate_fields(ϵ, a=rand(precision === :double ? Float64 : Float32, 20), b=2)
                
                # datatransfer.jl
                dev = dev_backend(dev_name)
                dev_mp2d, dev_mp3d, dev_grid2d, dev_grid3d, dev_attr = host2device(dev, mp2d, mp3d, grid2d, grid3d, attr)
                dev_mp2d_single = host2device(dev, mp2d)
                device2host!(mp2d, dev_mp2d)
                device2host!(mp2d, dev_mp2d, var=(:F,))
                device2host!(mp2d, dev_mp2d, var=(:F, :ξ, :σij))
                
                # calculator.jl - test deformation gradient updates
                T2_val = precision === :double ? 1.0 : 1.0f0
                dfs_2d = precision === :double ? [df1, df2, df3, df4] : Float32.([df1, df2, df3, df4])
                dfs_3d = precision === :double ? [df1, df2, df3, df4, df5, df6, df7, df8, df9] : Float32.([df1, df2, df3, df4, df5, df6, df7, df8, df9])
                
                update_F!(mp2d, dfs_2d[1], dfs_2d[2], dfs_2d[3], dfs_2d[4], T2_val, ix, Dim2D())
                update_F!(mp3d, dfs_3d[1], dfs_3d[2], dfs_3d[3], dfs_3d[4], dfs_3d[5], 
                         dfs_3d[6], dfs_3d[7], dfs_3d[8], dfs_3d[9], T2_val, ix, Dim3D())
                detF(mp2d, ix, Dim2D())
                detF(mp3d, ix, Dim3D())
                
                # material models - compile constitutive functions
                # Add material properties for constitutive model testing
                if hasfield(typeof(mp2d), :Ks)
                    material_mp2d = merge(mp2d, generate_fields(ϵ, Ks=rand(precision === :double ? Float64 : Float32, 5),
                        Gs=rand(precision === :double ? Float64 : Float32, 5),
                        c=rand(precision === :double ? Float64 : Float32, 5),
                        ϕ=rand(precision === :double ? Float64 : Float32, 5),
                        ψ=rand(precision === :double ? Float64 : Float32, 5),
                        σt=rand(precision === :double ? Float64 : Float32, 5)))
                    material_mp3d = merge(mp3d, generate_fields(ϵ, Ks=rand(precision === :double ? Float64 : Float32, 5),
                        Gs=rand(precision === :double ? Float64 : Float32, 5),
                        c=rand(precision === :double ? Float64 : Float32, 5),
                        ϕ=rand(precision === :double ? Float64 : Float32, 5),
                        ψ=rand(precision === :double ? Float64 : Float32, 5),
                        σt=rand(precision === :double ? Float64 : Float32, 5)))
                else
                    material_mp2d = merge(mp2d, (Ks=rand(precision === :double ? Float64 : Float32, 5),
                        Gs=rand(precision === :double ? Float64 : Float32, 5),
                        c=rand(precision === :double ? Float64 : Float32, 5),
                        ϕ=rand(precision === :double ? Float64 : Float32, 5),
                        ψ=rand(precision === :double ? Float64 : Float32, 5),
                        σt=rand(precision === :double ? Float64 : Float32, 5)))
                    material_mp3d = merge(mp3d, (Ks=rand(precision === :double ? Float64 : Float32, 5),
                        Gs=rand(precision === :double ? Float64 : Float32, 5),
                        c=rand(precision === :double ? Float64 : Float32, 5),
                        ϕ=rand(precision === :double ? Float64 : Float32, 5),
                        ψ=rand(precision === :double ? Float64 : Float32, 5),
                        σt=rand(precision === :double ? Float64 : Float32, 5)))
                end
                
                # interpolation.jl - test basis functions
                x1, x2 = precision === :double ? (0.5, 0.3) : (0.5f0, 0.3f0)
                invh = precision === :double ? 2.0 : 2.0f0
                linearbasis(x1, x2, invh)
            end
        end
    end
end