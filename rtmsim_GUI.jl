#include RTMsim
include("rtmsim.jl")
import .rtmsim

#include packages
using Gtk
using GLMakie
using Makie
using NativeFileDialog
using Gtk.ShortNames, Gtk.GConstants, Gtk.Graphics

#one Gtk window
win = GtkWindow("RTMsim"); 

#define buttons
sm=GtkButton("Select mesh file");pm=GtkButton("Plot mesh");ps=GtkButton("Plot sets")
ss=GtkButton("Start simulation");cs=GtkButton("Continue simulation")
sel=GtkButton("Select inlet port"); si=GtkButton("Start interactive");ci=GtkButton("Continue interactive");
sr=GtkButton("Select results file");pr=GtkButton("Plot results")
po=GtkButton("Plot overview")
pf=GtkButton("Plot filling")
q=GtkButton("Quit")
h=GtkButton("Help")
in1=GtkButton("Select input file")
in3=GtkButton("Run with input file")

#define input fields
in2=GtkEntry(); set_gtk_property!(in2,:text,"inputfiles\\input.txt");
mf=GtkEntry(); set_gtk_property!(mf,:text,"meshfiles\\mesh_permeameter1_foursets.bdf");
t=GtkEntry(); set_gtk_property!(t,:text,"200")
rf=GtkEntry(); set_gtk_property!(rf,:text,"results.jld2")
r=GtkEntry(); set_gtk_property!(r,:text,"0.01")
p1_0=GtkEntry(); set_gtk_property!(p1_0,:text,"Set 1");GAccessor.editable(GtkEditable(p1_0),false) 
p2_0=GtkEntry(); set_gtk_property!(p2_0,:text,"Set 2");GAccessor.editable(GtkEditable(p2_0),false) 
p3_0=GtkEntry(); set_gtk_property!(p3_0,:text,"Set 3");GAccessor.editable(GtkEditable(p3_0),false) 
p4_0=GtkEntry(); set_gtk_property!(p4_0,:text,"Set 4");GAccessor.editable(GtkEditable(p4_0),false) 
par_1=GtkEntry(); set_gtk_property!(par_1,:text,"135000")
par_2=GtkEntry(); set_gtk_property!(par_2,:text,"100000")
par_3=GtkEntry(); set_gtk_property!(par_3,:text,"0.06")
p0_1=GtkEntry(); set_gtk_property!(p0_1,:text,"0.003")
p0_2=GtkEntry(); set_gtk_property!(p0_2,:text,"0.7")
p0_3=GtkEntry(); set_gtk_property!(p0_3,:text,"3e-10")
p0_4=GtkEntry(); set_gtk_property!(p0_4,:text,"1")
p0_5=GtkEntry(); set_gtk_property!(p0_5,:text,"1")
p0_6=GtkEntry(); set_gtk_property!(p0_6,:text,"0")
p0_7=GtkEntry(); set_gtk_property!(p0_7,:text,"0")
p1_1=GtkEntry(); set_gtk_property!(p1_1,:text,"0.003")
p1_2=GtkEntry(); set_gtk_property!(p1_2,:text,"0.7")
p1_3=GtkEntry(); set_gtk_property!(p1_3,:text,"3e-10")
p1_4=GtkEntry(); set_gtk_property!(p1_4,:text,"1")
p1_5=GtkEntry(); set_gtk_property!(p1_5,:text,"1")
p1_6=GtkEntry(); set_gtk_property!(p1_6,:text,"0")
p1_7=GtkEntry(); set_gtk_property!(p1_7,:text,"0")
p2_1=GtkEntry(); set_gtk_property!(p2_1,:text,"0.003")
p2_2=GtkEntry(); set_gtk_property!(p2_2,:text,"0.7")
p2_3=GtkEntry(); set_gtk_property!(p2_3,:text,"3e-10")
p2_4=GtkEntry(); set_gtk_property!(p2_4,:text,"1")
p2_5=GtkEntry(); set_gtk_property!(p2_5,:text,"1")
p2_6=GtkEntry(); set_gtk_property!(p2_6,:text,"0")
p2_7=GtkEntry(); set_gtk_property!(p2_7,:text,"0")
p3_1=GtkEntry(); set_gtk_property!(p3_1,:text,"0.003")
p3_2=GtkEntry(); set_gtk_property!(p3_2,:text,"0.7")
p3_3=GtkEntry(); set_gtk_property!(p3_3,:text,"3e-10")
p3_4=GtkEntry(); set_gtk_property!(p3_4,:text,"1")
p3_5=GtkEntry(); set_gtk_property!(p3_5,:text,"1")
p3_6=GtkEntry(); set_gtk_property!(p3_6,:text,"0")
p3_7=GtkEntry(); set_gtk_property!(p3_7,:text,"0")
p4_1=GtkEntry(); set_gtk_property!(p4_1,:text,"0.003")
p4_2=GtkEntry(); set_gtk_property!(p4_2,:text,"0.7")
p4_3=GtkEntry(); set_gtk_property!(p4_3,:text,"3e-10")
p4_4=GtkEntry(); set_gtk_property!(p4_4,:text,"1")
p4_5=GtkEntry(); set_gtk_property!(p4_5,:text,"1")
p4_6=GtkEntry(); set_gtk_property!(p4_6,:text,"0")
p4_7=GtkEntry(); set_gtk_property!(p4_7,:text,"0")

#define radio buttons
choices = ["Ignore",  "Pressure inlet", "Pressure outlet", "Patch" ]
f1 = Gtk.GtkBox(:v);
r1 = Vector{RadioButton}(undef, 4)
r1[1] = RadioButton(choices[1]);                   push!(f1,r1[1])
r1[2] = RadioButton(r1[1],choices[2],active=true); push!(f1,r1[2])
r1[3] = RadioButton(r1[2],choices[3]);             push!(f1,r1[3])
r1[4] = RadioButton(r1[3],choices[4]);             push!(f1,r1[4])
f2 = Gtk.GtkBox(:v);
r2 = Vector{RadioButton}(undef, 4)
r2[1] = RadioButton(choices[1],active=true); push!(f2,r2[1])
r2[2] = RadioButton(r2[1],choices[2]);       push!(f2,r2[2])
r2[3] = RadioButton(r2[2],choices[3]);       push!(f2,r2[3])
r2[4] = RadioButton(r2[3],choices[4]);       push!(f2,r2[4])
f3 = Gtk.GtkBox(:v);
r3 = Vector{RadioButton}(undef, 4)
r3[1] = RadioButton(choices[1],active=true); push!(f3,r3[1])
r3[2] = RadioButton(r3[1],choices[2]);       push!(f3,r3[2])
r3[3] = RadioButton(r3[2],choices[3]);       push!(f3,r3[3])
r3[4] = RadioButton(r3[3],choices[4]);       push!(f3,r3[4])
f4 = Gtk.GtkBox(:v);
r4 = Vector{RadioButton}(undef, 4)
r4[1] = RadioButton(choices[1],active=true); push!(f4,r4[1])
r4[2] = RadioButton(r4[1],choices[2]);       push!(f4,r4[2])
r4[3] = RadioButton(r4[2],choices[3]);       push!(f4,r4[3])
r4[4] = RadioButton(r4[3],choices[4]);       push!(f4,r4[4])

#define logo
im=Gtk.GtkImage("figures\\rtmsim_logo1_h200px_grey.png")

#assembly elements in grid pattern
g = GtkGrid()    #Cartesian coordinates, g[column,row]
set_gtk_property!(g, :column_spacing, 5) 
set_gtk_property!(g, :row_spacing, 5) 
g[1,1]=sm; g[2,1] = mf; g[3,1] = pm;  g[4,1] = ps;  g[7,1] = in1;  g[8,1] = in2;  g[9,1] = in3;              
           g[2,2] = t;  g[3,2] = ss;  g[4,2] = cs; 
           g[2,3] = r;  g[3,3] = sel; g[4,3] = si; g[5,3] = ci; 
g[1,4] = sr; g[2,4] = rf; g[3,4] = pr; g[4,4] = po;g[5,4] = pf;
g[7:10,6:9] = im;
                                 g[3,11] = f1;   g[4,11] = f2;   g[5,11] = f3;   g[6,11] = f4;
g[1,12] = par_1; g[2,12] = p0_1; g[3,12] = p1_1; g[4,12] = p2_1; g[5,12] = p3_1; g[6,12] = p4_1; 
g[1,13] = par_2; g[2,13] = p0_2; g[3,13] = p1_2; g[4,13] = p2_2; g[5,13] = p3_2; g[6,13] = p4_2; 
g[1,14] = par_3; g[2,14] = p0_3; g[3,14] = p1_3; g[4,14] = p2_3; g[5,14] = p3_3; g[6,14] = p4_3; 
                 g[2,15] = p0_4; g[3,15] = p1_4; g[4,15] = p2_4; g[5,15] = p3_4; g[6,15] = p4_4; 
                 g[2,16] = p0_5; g[3,16] = p1_5; g[4,16] = p2_5; g[5,16] = p3_5; g[6,16] = p4_5; 
                 g[2,17] = p0_6; g[3,17] = p1_6; g[4,17] = p2_6; g[5,17] = p3_6; g[6,17] = p4_6;      g[9,17] = h; 
                 g[2,18] = p0_7; g[3,18] = p1_7; g[4,18] = p2_7; g[5,18] = p3_7; g[6,18] = p4_7;      g[9,18] = q; 
push!(win, g)

#callback functions
function sm_clicked(w)
    str = pick_file(pwd(),filterlist="bdf");
    set_gtk_property!(mf,:text,str);
end
function sr_clicked(w)
    str = pick_file(pwd(),filterlist="jld2");
    set_gtk_property!(rf,:text,str);
end
function pm_clicked(w)
    str = get_gtk_property(mf,:text,String)
    rtmsim.plot_mesh(str,1)
end
function ps_clicked(w)
    str = get_gtk_property(mf,:text,String)
    rtmsim.plot_sets(str)
end
function sel_clicked(w)
    str = get_gtk_property(mf,:text,String)
    rtmsim.plot_mesh(str,2)
end
function ss_clicked(w)
    str1 = get_gtk_property(mf,:text,String); str2 = get_gtk_property(t,:text,String); str3 = get_gtk_property(par_3,:text,String); str4 = get_gtk_property(par_1,:text,String); str5 = get_gtk_property(par_2,:text,String);
    str11 = get_gtk_property(p0_1,:text,String); str12 = get_gtk_property(p0_2,:text,String); str13 = get_gtk_property(p0_3,:text,String); str14 = get_gtk_property(p0_4,:text,String); str15 = get_gtk_property(p0_5,:text,String); str16 = get_gtk_property(p0_6,:text,String); str17 = get_gtk_property(p0_7,:text,String);
    str21 = get_gtk_property(p1_1,:text,String); str22 = get_gtk_property(p1_2,:text,String); str23 = get_gtk_property(p1_3,:text,String); str24 = get_gtk_property(p1_4,:text,String); str25 = get_gtk_property(p1_5,:text,String); str26 = get_gtk_property(p1_6,:text,String); str27 = get_gtk_property(p1_7,:text,String); 
    str31 = get_gtk_property(p2_1,:text,String); str32 = get_gtk_property(p2_2,:text,String); str33 = get_gtk_property(p2_3,:text,String); str34 = get_gtk_property(p2_4,:text,String); str35 = get_gtk_property(p2_5,:text,String); str36 = get_gtk_property(p2_6,:text,String); str37 = get_gtk_property(p2_7,:text,String);
    str41 = get_gtk_property(p3_1,:text,String); str42 = get_gtk_property(p3_2,:text,String); str43 = get_gtk_property(p3_3,:text,String); str44 = get_gtk_property(p3_4,:text,String); str45 = get_gtk_property(p3_5,:text,String); str46 = get_gtk_property(p3_6,:text,String); str47 = get_gtk_property(p3_7,:text,String); 
    str51 = get_gtk_property(p4_1,:text,String); str52 = get_gtk_property(p4_2,:text,String); str53 = get_gtk_property(p4_3,:text,String); str54 = get_gtk_property(p4_4,:text,String); str55 = get_gtk_property(p4_5,:text,String); str56 = get_gtk_property(p4_6,:text,String); str57 = get_gtk_property(p4_7,:text,String);
    if [get_gtk_property(b,:active,Bool) for b in r1] == [true, false, false, false]; patchtype1val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, true, false, false]; patchtype1val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, true, false]; patchtype1val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, false, true]; patchtype1val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r2] == [true, false, false, false]; patchtype2val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, true, false, false]; patchtype2val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, true, false]; patchtype2val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, false, true]; patchtype2val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r3] == [true, false, false, false]; patchtype3val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, true, false, false]; patchtype3val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, true, false]; patchtype3val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, false, true]; patchtype3val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r4] == [true, false, false, false]; patchtype4val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, true, false, false]; patchtype4val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, true, false]; patchtype4val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, false, true]; patchtype4val=Int64(2); end;
    str61="0.01"; #str61 = get_gtk_property(r,:text,String)
    restartval=Int64(0); interactiveval=Int64(0); noutval=Int64(16); 
    rtmsim.rtmsim_rev1(1,str1,parse(Float64,str2), 1.01325e5,1.225,1.4,parse(Float64,str3), parse(Float64,str4),parse(Float64,str5), parse(Float64,str11),
    parse(Float64,str12),parse(Float64,str13),parse(Float64,str14),parse(Float64,str15),parse(Float64,str16),parse(Float64,str17),
    parse(Float64,str21),parse(Float64,str22),parse(Float64,str23),parse(Float64,str24),parse(Float64,str25),parse(Float64,str26),parse(Float64,str27), 
    parse(Float64,str31),parse(Float64,str32),parse(Float64,str33),parse(Float64,str34),parse(Float64,str35),parse(Float64,str36),parse(Float64,str37),
    parse(Float64,str41),parse(Float64,str42),parse(Float64,str43),parse(Float64,str44),parse(Float64,str45),parse(Float64,str46),parse(Float64,str47),
    parse(Float64,str51),parse(Float64,str52),parse(Float64,str53),parse(Float64,str54),parse(Float64,str55),parse(Float64,str56),parse(Float64,str57),
    patchtype1val,patchtype2val,patchtype3val,patchtype4val, restartval,"results.jld2", interactiveval,parse(Float64,str61), noutval);
end
function cs_clicked(w)
    str1 = get_gtk_property(mf,:text,String); str2 = get_gtk_property(t,:text,String); str3 = get_gtk_property(par_3,:text,String); str4 = get_gtk_property(par_1,:text,String); str5 = get_gtk_property(par_2,:text,String);
    str11 = get_gtk_property(p0_1,:text,String); str12 = get_gtk_property(p0_2,:text,String); str13 = get_gtk_property(p0_3,:text,String); str14 = get_gtk_property(p0_4,:text,String); str15 = get_gtk_property(p0_5,:text,String); str16 = get_gtk_property(p0_6,:text,String); str17 = get_gtk_property(p0_7,:text,String);
    str21 = get_gtk_property(p1_1,:text,String); str22 = get_gtk_property(p1_2,:text,String); str23 = get_gtk_property(p1_3,:text,String); str24 = get_gtk_property(p1_4,:text,String); str25 = get_gtk_property(p1_5,:text,String); str26 = get_gtk_property(p1_6,:text,String); str27 = get_gtk_property(p1_7,:text,String); 
    str31 = get_gtk_property(p2_1,:text,String); str32 = get_gtk_property(p2_2,:text,String); str33 = get_gtk_property(p2_3,:text,String); str34 = get_gtk_property(p2_4,:text,String); str35 = get_gtk_property(p2_5,:text,String); str36 = get_gtk_property(p2_6,:text,String); str37 = get_gtk_property(p2_7,:text,String);
    str41 = get_gtk_property(p3_1,:text,String); str42 = get_gtk_property(p3_2,:text,String); str43 = get_gtk_property(p3_3,:text,String); str44 = get_gtk_property(p3_4,:text,String); str45 = get_gtk_property(p3_5,:text,String); str46 = get_gtk_property(p3_6,:text,String); str47 = get_gtk_property(p3_7,:text,String); 
    str51 = get_gtk_property(p4_1,:text,String); str52 = get_gtk_property(p4_2,:text,String); str53 = get_gtk_property(p4_3,:text,String); str54 = get_gtk_property(p4_4,:text,String); str55 = get_gtk_property(p4_5,:text,String); str56 = get_gtk_property(p4_6,:text,String); str57 = get_gtk_property(p4_7,:text,String);
    if [get_gtk_property(b,:active,Bool) for b in r1] == [true, false, false, false]; patchtype1val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, true, false, false]; patchtype1val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, true, false]; patchtype1val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, false, true]; patchtype1val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r2] == [true, false, false, false]; patchtype2val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, true, false, false]; patchtype2val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, true, false]; patchtype2val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, false, true]; patchtype2val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r3] == [true, false, false, false]; patchtype3val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, true, false, false]; patchtype3val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, true, false]; patchtype3val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, false, true]; patchtype3val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r4] == [true, false, false, false]; patchtype4val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, true, false, false]; patchtype4val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, true, false]; patchtype4val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, false, true]; patchtype4val=Int64(2); end;
    str61="0.01"; #str61 = get_gtk_property(r,:text,String)
    restartval=Int64(1); interactiveval=Int64(0); noutval=Int64(16); 
    rtmsim.rtmsim_rev1(1,str1,parse(Float64,str2), 1.01325e5,1.225,1.4,parse(Float64,str3), parse(Float64,str4),parse(Float64,str5), parse(Float64,str11),
    parse(Float64,str12),parse(Float64,str13),parse(Float64,str14),parse(Float64,str15),parse(Float64,str16),parse(Float64,str17),
    parse(Float64,str21),parse(Float64,str22),parse(Float64,str23),parse(Float64,str24),parse(Float64,str25),parse(Float64,str26),parse(Float64,str27), 
    parse(Float64,str31),parse(Float64,str32),parse(Float64,str33),parse(Float64,str34),parse(Float64,str35),parse(Float64,str36),parse(Float64,str37),
    parse(Float64,str41),parse(Float64,str42),parse(Float64,str43),parse(Float64,str44),parse(Float64,str45),parse(Float64,str46),parse(Float64,str47),
    parse(Float64,str51),parse(Float64,str52),parse(Float64,str53),parse(Float64,str54),parse(Float64,str55),parse(Float64,str56),parse(Float64,str57),
    patchtype1val,patchtype2val,patchtype3val,patchtype4val, restartval,"results.jld2", interactiveval,parse(Float64,str61), noutval);
end
function si_clicked(w)
    str1 = get_gtk_property(mf,:text,String); str2 = get_gtk_property(t,:text,String); str3 = get_gtk_property(par_3,:text,String); str4 = get_gtk_property(par_1,:text,String); str5 = get_gtk_property(par_2,:text,String);
    str11 = get_gtk_property(p0_1,:text,String); str12 = get_gtk_property(p0_2,:text,String); str13 = get_gtk_property(p0_3,:text,String); str14 = get_gtk_property(p0_4,:text,String); str15 = get_gtk_property(p0_5,:text,String); str16 = get_gtk_property(p0_6,:text,String); str17 = get_gtk_property(p0_7,:text,String);
    str21 = get_gtk_property(p1_1,:text,String); str22 = get_gtk_property(p1_2,:text,String); str23 = get_gtk_property(p1_3,:text,String); str24 = get_gtk_property(p1_4,:text,String); str25 = get_gtk_property(p1_5,:text,String); str26 = get_gtk_property(p1_6,:text,String); str27 = get_gtk_property(p1_7,:text,String); 
    str31 = get_gtk_property(p2_1,:text,String); str32 = get_gtk_property(p2_2,:text,String); str33 = get_gtk_property(p2_3,:text,String); str34 = get_gtk_property(p2_4,:text,String); str35 = get_gtk_property(p2_5,:text,String); str36 = get_gtk_property(p2_6,:text,String); str37 = get_gtk_property(p2_7,:text,String);
    str41 = get_gtk_property(p3_1,:text,String); str42 = get_gtk_property(p3_2,:text,String); str43 = get_gtk_property(p3_3,:text,String); str44 = get_gtk_property(p3_4,:text,String); str45 = get_gtk_property(p3_5,:text,String); str46 = get_gtk_property(p3_6,:text,String); str47 = get_gtk_property(p3_7,:text,String); 
    str51 = get_gtk_property(p4_1,:text,String); str52 = get_gtk_property(p4_2,:text,String); str53 = get_gtk_property(p4_3,:text,String); str54 = get_gtk_property(p4_4,:text,String); str55 = get_gtk_property(p4_5,:text,String); str56 = get_gtk_property(p4_6,:text,String); str57 = get_gtk_property(p4_7,:text,String);
    if [get_gtk_property(b,:active,Bool) for b in r1] == [true, false, false, false]; patchtype1val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, true, false, false]; patchtype1val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, true, false]; patchtype1val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, false, true]; patchtype1val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r2] == [true, false, false, false]; patchtype2val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, true, false, false]; patchtype2val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, true, false]; patchtype2val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, false, true]; patchtype2val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r3] == [true, false, false, false]; patchtype3val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, true, false, false]; patchtype3val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, true, false]; patchtype3val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, false, true]; patchtype3val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r4] == [true, false, false, false]; patchtype4val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, true, false, false]; patchtype4val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, true, false]; patchtype4val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, false, true]; patchtype4val=Int64(2); end;
    str61 = get_gtk_property(r,:text,String)
    restartval=Int64(0); interactiveval=Int64(2); noutval=Int64(16); 
    rtmsim.rtmsim_rev1(1,str1,parse(Float64,str2), 1.01325e5,1.225,1.4,parse(Float64,str3), parse(Float64,str4),parse(Float64,str5), parse(Float64,str11),
    parse(Float64,str12),parse(Float64,str13),parse(Float64,str14),parse(Float64,str15),parse(Float64,str16),parse(Float64,str17),
    parse(Float64,str21),parse(Float64,str22),parse(Float64,str23),parse(Float64,str24),parse(Float64,str25),parse(Float64,str26),parse(Float64,str27), 
    parse(Float64,str31),parse(Float64,str32),parse(Float64,str33),parse(Float64,str34),parse(Float64,str35),parse(Float64,str36),parse(Float64,str37),
    parse(Float64,str41),parse(Float64,str42),parse(Float64,str43),parse(Float64,str44),parse(Float64,str45),parse(Float64,str46),parse(Float64,str47),
    parse(Float64,str51),parse(Float64,str52),parse(Float64,str53),parse(Float64,str54),parse(Float64,str55),parse(Float64,str56),parse(Float64,str57),
    patchtype1val,patchtype2val,patchtype3val,patchtype4val, restartval,"results.jld2", interactiveval,parse(Float64,str61), noutval);
end
function ci_clicked(w)
    str1 = get_gtk_property(mf,:text,String); str2 = get_gtk_property(t,:text,String); str3 = get_gtk_property(par_3,:text,String); str4 = get_gtk_property(par_1,:text,String); str5 = get_gtk_property(par_2,:text,String);
    str11 = get_gtk_property(p0_1,:text,String); str12 = get_gtk_property(p0_2,:text,String); str13 = get_gtk_property(p0_3,:text,String); str14 = get_gtk_property(p0_4,:text,String); str15 = get_gtk_property(p0_5,:text,String); str16 = get_gtk_property(p0_6,:text,String); str17 = get_gtk_property(p0_7,:text,String);
    str21 = get_gtk_property(p1_1,:text,String); str22 = get_gtk_property(p1_2,:text,String); str23 = get_gtk_property(p1_3,:text,String); str24 = get_gtk_property(p1_4,:text,String); str25 = get_gtk_property(p1_5,:text,String); str26 = get_gtk_property(p1_6,:text,String); str27 = get_gtk_property(p1_7,:text,String); 
    str31 = get_gtk_property(p2_1,:text,String); str32 = get_gtk_property(p2_2,:text,String); str33 = get_gtk_property(p2_3,:text,String); str34 = get_gtk_property(p2_4,:text,String); str35 = get_gtk_property(p2_5,:text,String); str36 = get_gtk_property(p2_6,:text,String); str37 = get_gtk_property(p2_7,:text,String);
    str41 = get_gtk_property(p3_1,:text,String); str42 = get_gtk_property(p3_2,:text,String); str43 = get_gtk_property(p3_3,:text,String); str44 = get_gtk_property(p3_4,:text,String); str45 = get_gtk_property(p3_5,:text,String); str46 = get_gtk_property(p3_6,:text,String); str47 = get_gtk_property(p3_7,:text,String); 
    str51 = get_gtk_property(p4_1,:text,String); str52 = get_gtk_property(p4_2,:text,String); str53 = get_gtk_property(p4_3,:text,String); str54 = get_gtk_property(p4_4,:text,String); str55 = get_gtk_property(p4_5,:text,String); str56 = get_gtk_property(p4_6,:text,String); str57 = get_gtk_property(p4_7,:text,String);
    if [get_gtk_property(b,:active,Bool) for b in r1] == [true, false, false, false]; patchtype1val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, true, false, false]; patchtype1val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, true, false]; patchtype1val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r1] == [false, false, false, true]; patchtype1val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r2] == [true, false, false, false]; patchtype2val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, true, false, false]; patchtype2val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, true, false]; patchtype2val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r2] == [false, false, false, true]; patchtype2val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r3] == [true, false, false, false]; patchtype3val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, true, false, false]; patchtype3val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, true, false]; patchtype3val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r3] == [false, false, false, true]; patchtype3val=Int64(2); end;
    if [get_gtk_property(b,:active,Bool) for b in r4] == [true, false, false, false]; patchtype4val=Int64(0);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, true, false, false]; patchtype4val=Int64(1);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, true, false]; patchtype4val=Int64(3);    elseif [get_gtk_property(b,:active,Bool) for b in r4] == [false, false, false, true]; patchtype4val=Int64(2); end;
    str61 = get_gtk_property(r,:text,String)
    restartval=Int64(1); interactiveval=Int64(2); noutval=Int64(16); 
    rtmsim.rtmsim_rev1(1,str1,parse(Float64,str2), 1.01325e5,1.225,1.4,parse(Float64,str3), parse(Float64,str4),parse(Float64,str5), parse(Float64,str11),
    parse(Float64,str12),parse(Float64,str13),parse(Float64,str14),parse(Float64,str15),parse(Float64,str16),parse(Float64,str17),
    parse(Float64,str21),parse(Float64,str22),parse(Float64,str23),parse(Float64,str24),parse(Float64,str25),parse(Float64,str26),parse(Float64,str27), 
    parse(Float64,str31),parse(Float64,str32),parse(Float64,str33),parse(Float64,str34),parse(Float64,str35),parse(Float64,str36),parse(Float64,str37),
    parse(Float64,str41),parse(Float64,str42),parse(Float64,str43),parse(Float64,str44),parse(Float64,str45),parse(Float64,str46),parse(Float64,str47),
    parse(Float64,str51),parse(Float64,str52),parse(Float64,str53),parse(Float64,str54),parse(Float64,str55),parse(Float64,str56),parse(Float64,str57),
    patchtype1val,patchtype2val,patchtype3val,patchtype4val, restartval,"results.jld2", interactiveval,parse(Float64,str61), noutval);
end
function pr_clicked(w)
    str = get_gtk_property(rf,:text,String)
    rtmsim.plot_results(str) 
end
function po_clicked(w)
    rtmsim.plot_overview(-1,16) 
end
function pf_clicked(w)
    rtmsim.plot_filling(-1,16) 
end
function q_clicked(w)
    GLMakie.destroy!(GLMakie.global_gl_screen())
    Gtk.destroy(win)
end
function h_clicked(w)
    i=GtkImage("figures\\rtmsim_help.png");
    w=GtkWindow(i,"Help");
    show(i);
end
function in1_clicked(w)
    str = pick_file(pwd(),filterlist="txt");
    set_gtk_property!(in2,:text,str);
end
function in3_clicked(w)
    str = get_gtk_property(in2,:text,String)
    rtmsim.start_rtmsim(str)
end

#callbacks
signal_connect(sm_clicked,sm,"clicked")
signal_connect(sr_clicked,sr,"clicked")
signal_connect(pm_clicked,pm,"clicked")
signal_connect(ps_clicked,ps,"clicked")
signal_connect(ss_clicked,ss,"clicked");
signal_connect(cs_clicked,cs,"clicked")
signal_connect(sel_clicked,sel,"clicked")
signal_connect(si_clicked,si,"clicked");
signal_connect(ci_clicked,ci,"clicked")
signal_connect(pr_clicked,pr,"clicked")
signal_connect(po_clicked,po,"clicked")
signal_connect(pf_clicked,pf,"clicked")
signal_connect(q_clicked,q,"clicked")
signal_connect(h_clicked,h,"clicked")
signal_connect(in1_clicked,in1,"clicked")
signal_connect(in3_clicked,in3,"clicked")

#show GUI
showall(win);
if !isinteractive()
    c = Condition()
    signal_connect(win, :destroy) do widget
        notify(c)
    end
    wait(c)
end
