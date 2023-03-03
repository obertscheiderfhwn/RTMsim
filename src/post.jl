"""
    function plot_results(resultsfilename)

Create contour plots of the filling factor and the pressure after loading a results file.

Arguments:
- resultsfilename :: String

Unit test: 
- `WORK_DIR=pwd(); resultsfilename=joinpath(WORK_DIR,"results.jld2"); rtmsim.plot_results(resultsfilename);`
"""
function plot_results(resultsfilename)
    #create contour plots of the filling factor and the pressure after loading a results file
    #default call: rtmsim.plot_results("results.jld2")

    if ~isfile(resultsfilename)
        errorstring=string("File ",resultsfilename," not existing"* "\n") 
        error(errorstring)
    end
    t_digits=2 
    t_div=10^2
    @load resultsfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out

    gamma_plot=Vector{Float64}(undef, N)
    deltap=maximum(p_new)-minimum(p_new)     
    for ind=1:N
        if gamma_out[ind]>0.8
            gamma_plot[ind]=1
        else
            gamma_plot[ind]=0
        end
    end
    deltagamma=1  #deltagamma=maximum(gamma_plot)-minimum(gamma_plot)

    #for poly plot
    inds0=findall(gamma_out.>-0.5)
    N0=length(inds0)
    X=Array{Float64}(undef, 3, N0)
    Y=Array{Float64}(undef, 3, N0)
    Z=Array{Float64}(undef, 3, N0)
    C_p=Array{Float32}(undef, 3, N0)        
    C_gamma=Array{Float32}(undef, 3, N0)
    inds1=findall(gamma_out.<-0.5)
    N1=length(inds1)
    X1=Array{Float64}(undef, 3, N1)
    Y1=Array{Float64}(undef, 3, N1)
    Z1=Array{Float64}(undef, 3, N1)  
    C1_gamma=Array{Float32}(undef, 3, N1)
    C1_p=Array{Float32}(undef, 3, N1)
    for i in 1:N0
        ind=inds0[i]
        X[1,i]=gridx[cellgridid[ind,1]]
        X[2,i]=gridx[cellgridid[ind,2]]
        X[3,i]=gridx[cellgridid[ind,3]]
        Y[1,i]=gridy[cellgridid[ind,1]]
        Y[2,i]=gridy[cellgridid[ind,2]]
        Y[3,i]=gridy[cellgridid[ind,3]]
        Z[1,i]=gridz[cellgridid[ind,1]]
        Z[2,i]=gridz[cellgridid[ind,2]]
        Z[3,i]=gridz[cellgridid[ind,3]]
        C_gamma[1,i]=gamma_plot[ind]/deltagamma
        C_gamma[2,i]=gamma_plot[ind]/deltagamma
        C_gamma[3,i]=gamma_plot[ind]/deltagamma
        C_p[1,i]=p_new[ind]/deltap
        C_p[2,i]=p_new[ind]/deltap
        C_p[3,i]=p_new[ind]/deltap
    end
    xyz = reshape([X[:] Y[:] Z[:]]', :)        
    for i in 1:N1
        ind=inds1[i]
        X1[1,i]=gridx[cellgridid[ind,1]]
        X1[2,i]=gridx[cellgridid[ind,2]]
        X1[3,i]=gridx[cellgridid[ind,3]]
        Y1[1,i]=gridy[cellgridid[ind,1]]
        Y1[2,i]=gridy[cellgridid[ind,2]]
        Y1[3,i]=gridy[cellgridid[ind,3]]
        Z1[1,i]=gridz[cellgridid[ind,1]]
        Z1[2,i]=gridz[cellgridid[ind,2]]
        Z1[3,i]=gridz[cellgridid[ind,3]]
        C1_gamma[1,i]=0.5
        C1_gamma[2,i]=0.5
        C1_gamma[3,i]=0.5
        C1_p[1,i]=p_new[ind]/deltap
        C1_p[2,i]=p_new[ind]/deltap
        C1_p[3,i]=p_new[ind]/deltap
    end
    xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

    #bounding box
    deltax=maximum(gridx)-minimum(gridx)
    deltay=maximum(gridy)-minimum(gridy)
    deltaz=maximum(gridz)-minimum(gridz)
    mindelta=min(deltax,deltay,deltaz)
    maxdelta=max(deltax,deltay,deltaz)
    if mindelta<maxdelta*0.001
        eps_delta=maxdelta*0.001
    else
        eps_delta=0
    end 
    ax=(deltax+eps_delta)/(mindelta+eps_delta)
    ay=(deltay+eps_delta)/(mindelta+eps_delta)
    az=(deltaz+eps_delta)/(mindelta+eps_delta)

    resolution_val=600
    fig = Figure(resolution=(2*resolution_val, resolution_val))
    ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
    poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
    if N1>0 
        poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
    end
    hidedecorations!(ax1)
    hidespines!(ax1) 
    ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Pressure at t=", string(round(t_div*t)/t_div) ,"s"))
    poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_p[:], strokewidth=1, colorrange=(0,1))
    if N1>0 
        poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_p[:], strokewidth=1, colorrange=(0,1))
    end
    Colorbar(fig[1, 3], limits = (0, deltap), colormap = :viridis,  vertical=true, height=Relative(0.5))  
    #Colorbar(fig[2, 2], limits = (0, deltap), colormap = :viridis,  vertical=false, width=Relative(0.5))
    hidedecorations!(ax2)
    hidespines!(ax2) 
    display(fig)
end 


"""
    function plot_overview(n_out,n_pics)

Create filling factor contour plots. `n_out` is the index of the last output file, if `n_out==-1` the output file with the highest index is chosen. Consider the last `n_pics` for creating the contour plots at four equidistant time intervals, if `n_pics==-1` all available output files are considered.

Arguments:
- n_out :: Int
- n_pics :: Int

Unit test: 
- `rtmsim.plot_overview(-1,-1)`
"""
function plot_overview(n_out,n_pics)
    #

    val=0
    n_out_start=-1
    if n_out==-1
        vec1=glob("output_*.jld2")
        for i=1:length(vec1)
            vec2=split(vec1[i],".")
            vec3=split(vec2[1],"_")
            val=max(val,parse(Int64,vec3[2]))
            if i==1
                n_out_start=parse(Int64,vec3[2])
            end
        end            
        n_out=val
    end
    if n_pics==-1
        n_pics=(n_out-n_out_start)
    end
    if mod(n_pics,4)!=0
        errorstring=string("n_pics must be multiple of four"* "\n") 
        error(errorstring)
    end
    t_digits=2 
    t_div=10^2

    resolution_val=300
    fig = Figure(resolution=(4*resolution_val, resolution_val))
    i_out=n_out-3*Int64(n_pics/4)
    for i_plot in 1:4          
        outputfilename=string("output_",string(i_out),".jld2")
        if ~isfile(outputfilename)
            errorstring=string("File ",outputfilename," not existing"* "\n") 
            error(errorstring)
        else
            loadfilename="results_temp.jld2"
            cp(outputfilename,loadfilename;force=true)
            @load loadfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
        end

        gamma_plot=Vector{Float64}(undef, N)
        deltap=maximum(p_new)-minimum(p_new)      
        for ind=1:N
            if gamma_out[ind]>0.8
                gamma_plot[ind]=1
            else
                gamma_plot[ind]=0
            end
        end
        deltagamma=1  #deltagamma=maximum(gamma_plot)-minimum(gamma_plot)

        #for poly plot
        inds0=findall(gamma_out.>-0.5)
        N0=length(inds0)
        X=Array{Float64}(undef, 3, N0)
        Y=Array{Float64}(undef, 3, N0)
        Z=Array{Float64}(undef, 3, N0)
        C_p=Array{Float32}(undef, 3, N0)        
        C_gamma=Array{Float32}(undef, 3, N0)
        inds1=findall(gamma_out.<-0.5)
        N1=length(inds1)
        X1=Array{Float64}(undef, 3, N1)
        Y1=Array{Float64}(undef, 3, N1)
        Z1=Array{Float64}(undef, 3, N1)  
        C1_p=Array{Float32}(undef, 3, N1)
        C1_gamma=Array{Float32}(undef, 3, N1)
        for i in 1:N0
            ind=inds0[i]
            X[1,i]=gridx[cellgridid[ind,1]]
            X[2,i]=gridx[cellgridid[ind,2]]
            X[3,i]=gridx[cellgridid[ind,3]]
            Y[1,i]=gridy[cellgridid[ind,1]]
            Y[2,i]=gridy[cellgridid[ind,2]]
            Y[3,i]=gridy[cellgridid[ind,3]]
            Z[1,i]=gridz[cellgridid[ind,1]]
            Z[2,i]=gridz[cellgridid[ind,2]]
            Z[3,i]=gridz[cellgridid[ind,3]]
            C_gamma[1,i]=gamma_plot[ind]/deltagamma
            C_gamma[2,i]=gamma_plot[ind]/deltagamma
            C_gamma[3,i]=gamma_plot[ind]/deltagamma
            C_p[1,i]=p_new[ind]/deltap
            C_p[2,i]=p_new[ind]/deltap
            C_p[3,i]=p_new[ind]/deltap
        end
        xyz = reshape([X[:] Y[:] Z[:]]', :)
        for i in 1:N1
            ind=inds1[i]
            X1[1,i]=gridx[cellgridid[ind,1]]
            X1[2,i]=gridx[cellgridid[ind,2]]
            X1[3,i]=gridx[cellgridid[ind,3]]
            Y1[1,i]=gridy[cellgridid[ind,1]]
            Y1[2,i]=gridy[cellgridid[ind,2]]
            Y1[3,i]=gridy[cellgridid[ind,3]]
            Z1[1,i]=gridz[cellgridid[ind,1]]
            Z1[2,i]=gridz[cellgridid[ind,2]]
            Z1[3,i]=gridz[cellgridid[ind,3]]
            C1_gamma[1,i]=0.5
            C1_gamma[2,i]=0.5
            C1_gamma[3,i]=0.5
            C1_p[1,i]=p_new[ind]/deltap
            C1_p[2,i]=p_new[ind]/deltap
            C1_p[3,i]=p_new[ind]/deltap
        end
        xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

        #bounding box
        deltax=maximum(gridx)-minimum(gridx)
        deltay=maximum(gridy)-minimum(gridy)
        deltaz=maximum(gridz)-minimum(gridz)
        mindelta=min(deltax,deltay,deltaz)
        maxdelta=max(deltax,deltay,deltaz)
        if mindelta<maxdelta*0.001
            eps_delta=maxdelta*0.001
        else
            eps_delta=0
        end 
        ax=(deltax+eps_delta)/(mindelta+eps_delta)
        ay=(deltay+eps_delta)/(mindelta+eps_delta)
        az=(deltaz+eps_delta)/(mindelta+eps_delta)
        if i_plot==1
            ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
            poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
            if N1>0
                poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
            end
            hidedecorations!(ax1)
            hidespines!(ax1) 
        elseif i_plot==2
            ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
            poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
            if N1>0
                poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
            end
            hidedecorations!(ax2)
            hidespines!(ax2) 
        elseif i_plot==3
            ax3 = Axis3(fig[1, 3]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
            poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
            if N1>0
                poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
            end
            hidedecorations!(ax3)
            hidespines!(ax3) 
        elseif i_plot==4
            ax4 = Axis3(fig[1, 4]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
            poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
            if N1>0
                poly!(connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
            end
            hidedecorations!(ax4)
            hidespines!(ax4) 
        end
        i_out=i_out+Int64(n_pics/4)
    end
    display(fig)
end 


"""
    function plot_filling(n_out,n_pics)

Create a window showing the filling factor contour plot at a selected time instance. Selection is with slider bar. `n_out` is the index of the last output file, if `n_out==-1` the output file with the highest index is chosen. Consider the last `n_pics` for creating the contour plots at four equidistant time intervals, if `n_pics==-1` all available output files are considered.

Arguments: 
- n_out :: Int
- n_pics :: Int

Unit test:
- `rtmsim.plot_filling(-1,-1)` 
"""
function plot_filling(n_out,n_pics)
    val=0
    n_out_start=-1
    if n_out==-1
        vec1=glob("output_*.jld2")
        for i=1:length(vec1)
            vec2=split(vec1[i],".")
            vec3=split(vec2[1],"_")
            val=max(val,parse(Int64,vec3[2]))
            if i==1
                n_out_start=parse(Int64,vec3[2])
            end
        end            
        n_out=val
    end
    if n_pics==-1
        n_pics=(n_out-n_out_start)
    end
    #plot the last n_pics pictures
    if n_pics<4
        errorstring=string("Makes no sense for n_pics<4"* "\n") 
        error(errorstring)
    end
    t_digits=2 
    t_div=10^2

    time_vector=[]
    output_array=[]
    inds=[]
    inds0=[]
    inds1=[]
    N=Int64(0)
    N0=Int64(0)
    N1=Int64(0)
    ax=0.0
    ay=0.0
    az=0.0
    t=0.0
    xyz=Vector{Float64}
    xyz1=Vector{Float64}
    X=Vector{Float64}
    Y=Vector{Float64}
    Z=Vector{Float64}
    C=Vector{Float64}
    X1=Vector{Float64}
    Y1=Vector{Float64}
    Z1=Vector{Float64}
    C1_gamma=Vector{Float64}        
    i_out=n_out-n_pics
    i_firstfile=1
    for i_plot in 1:n_pics+1          
        outputfilename=string("output_",string(i_out),".jld2")
        if ~isfile(outputfilename)
            errorstring=string("File ",outputfilename," not existing"* "\n") 
            error(errorstring)
        else
            loadfilename="results_temp.jld2"
            cp(outputfilename,loadfilename;force=true)
            @load loadfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
            if i_firstfile==1
                i_firstfile=0
                #for poly plot
                inds0=findall(gamma_out.>-0.5)
                N0=length(inds0)
                X=Array{Float64}(undef, 3, N0)
                Y=Array{Float64}(undef, 3, N0)
                Z=Array{Float64}(undef, 3, N0)
                C_p=Array{Float32}(undef, 3, N0)        
                C_gamma=Array{Float32}(undef, 3, N0)
                inds1=findall(gamma_out.<-0.5)
                N1=length(inds1)
                X1=Array{Float64}(undef, 3, N1)
                Y1=Array{Float64}(undef, 3, N1)
                Z1=Array{Float64}(undef, 3, N1)  
                C1_gamma=Array{Float32}(undef, 3, N1)
                for i in 1:N0
                    ind=inds0[i]
                    X[1,i]=gridx[cellgridid[ind,1]]
                    X[2,i]=gridx[cellgridid[ind,2]]
                    X[3,i]=gridx[cellgridid[ind,3]]
                    Y[1,i]=gridy[cellgridid[ind,1]]
                    Y[2,i]=gridy[cellgridid[ind,2]]
                    Y[3,i]=gridy[cellgridid[ind,3]]
                    Z[1,i]=gridz[cellgridid[ind,1]]
                    Z[2,i]=gridz[cellgridid[ind,2]]
                    Z[3,i]=gridz[cellgridid[ind,3]]
                end
                xyz = reshape([X[:] Y[:] Z[:]]', :)
                for i in 1:N1
                    ind=inds1[i]
                    X1[1,i]=gridx[cellgridid[ind,1]]
                    X1[2,i]=gridx[cellgridid[ind,2]]
                    X1[3,i]=gridx[cellgridid[ind,3]]
                    Y1[1,i]=gridy[cellgridid[ind,1]]
                    Y1[2,i]=gridy[cellgridid[ind,2]]
                    Y1[3,i]=gridy[cellgridid[ind,3]]
                    Z1[1,i]=gridz[cellgridid[ind,1]]
                    Z1[2,i]=gridz[cellgridid[ind,2]]
                    Z1[3,i]=gridz[cellgridid[ind,3]]
                    C1_gamma[1,i]=0.5
                    C1_gamma[2,i]=0.5
                    C1_gamma[3,i]=0.5
                end
                xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

                #bounding box
                deltax=maximum(gridx)-minimum(gridx)
                deltay=maximum(gridy)-minimum(gridy)
                deltaz=maximum(gridz)-minimum(gridz)
                mindelta=min(deltax,deltay,deltaz)
                maxdelta=max(deltax,deltay,deltaz)
                if mindelta<maxdelta*0.001
                    eps_delta=maxdelta*0.001
                else
                    eps_delta=0
                end 
                ax=(deltax+eps_delta)/(mindelta+eps_delta)
                ay=(deltay+eps_delta)/(mindelta+eps_delta)
                az=(deltaz+eps_delta)/(mindelta+eps_delta)
                time_vector=t
                output_array=gamma_out 
                N_val=N

            else
                time_vector=vcat(time_vector,t)
                output_array=hcat(output_array,gamma_out)
            end
        end
        i_out=i_out+1
    end

    gamma_plot=output_array[:,end]  
    for ind=1:N
        if gamma_plot[ind]>0.8
            gamma_plot[ind]=1
        else
            gamma_plot[ind]=0
        end
    end
    deltagamma=1  #deltagamma=maximum(gamma_plot)-minimum(gamma_plot)

    C_gamma=Array{Float32}(undef, 3, N0)
    for i in 1:N0
        ind=inds0[i]
        C_gamma[1,i]=gamma_plot[ind]/deltagamma
        C_gamma[2,i]=gamma_plot[ind]/deltagamma
        C_gamma[3,i]=gamma_plot[ind]/deltagamma
    end

    resolution_val=600
    fig = Figure(resolution=(resolution_val, resolution_val))   
    ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
    #p1=poly!(ax1,connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
    #if N1>0
    #    p2=poly!(ax1,connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
    #end
    hidedecorations!(ax1)
    hidespines!(ax1) 
    #display(fig)
    #sl_t = Slider(fig[2, 1], range = time_vector[1]:  (time_vector[end]-time_vector[1])/n_pics :time_vector[end], startvalue =  time_vector[end] )
    sl_t = Slider(fig[2, 1], range = time_vector[1]:  (time_vector[end]-time_vector[1])/n_pics :time_vector[end], startvalue =  time_vector[1] )
    point = lift(sl_t.value) do x           
        if x<0.5*(time_vector[end]+time_vector[1])
            gamma_plot=output_array[:,1]
        else
            gamma_plot=output_array[:,end]
        end
        time_val=x
        tind=1
        tind1=1
        tind2=2
        for i in 1:length(time_vector)-1
            if x>=0.5*(time_vector[i]+time_vector[i+1])
                tind=i+1
            end
        end            
        gamma_plot=output_array[:,tind]
        time_val=time_vector[tind]

        for ind=1:N
            if gamma_plot[ind]>0.8
                gamma_plot[ind]=1
            else
                gamma_plot[ind]=0
            end
        end
        deltagamma=1  #maximum(gamma_plot)-minimum(gamma_plot)
        for i in 1:N0
            ind=inds0[i]
            C_gamma[1,i]=gamma_plot[ind]/deltagamma
            C_gamma[2,i]=gamma_plot[ind]/deltagamma
            C_gamma[3,i]=gamma_plot[ind]/deltagamma
        end
        empty!(ax1.scene)
        p1=poly!(ax1,connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
        if N1>0
            p2=poly!(ax1,connect(xyz1, Makie.Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
        end
        ax1.title=string("Filling factor at t=", string(round(t_div*time_val)/t_div) ,"s")
        hidedecorations!(ax1)
        hidespines!(ax1) 
        display(fig)
    end
    
end 
