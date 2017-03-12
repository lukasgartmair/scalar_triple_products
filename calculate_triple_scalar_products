#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 22:35:34 2017

@author: lukas
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#https://en.wikipedia.org/wiki/Parallelepiped
#If a, b, and c are the parallelepiped edge lengths, and α, β, and γ are the internal angles between the edges, the volume is
#    V = a b c 1 + 2 cos ⁡ ( α ) cos ⁡ ( β ) cos ⁡ ( γ ) − cos 2 ⁡ ( α ) − cos 2 ⁡ ( β ) − cos 2 ⁡ ( γ ) .
# {\displaystyle V=abc{\sqrt {1+2\cos(\alpha )\cos(\beta )\cos(\gamma )-\cos ^{2}(\alpha )-\cos ^{2}(\beta )-\cos ^{2}(\gamma )\,}}.} 
# V=abc{\sqrt {1+2\cos(\alpha )\cos(\beta )\cos(\gamma )-\cos ^{2}(\alpha )-\cos ^{2}(\beta )-\cos ^{2}(\gamma )\,}}. 



lineStyles = {':': '_draw_dotted', '--': '_draw_dashed', '-.': '_draw_dash_dot', '-': '_draw_solid'}
line_styles_list = list(lineStyles.keys())


plot_types = ['sharp', 'dist', 'real']
filenames = ['sharp_ref', 'dist_c60', 'real']

for i,p in enumerate(plot_types):

        
        vdb_verts = np.genfromtxt(p + '_vdb_vertices.txt')
        vdb_faces = np.genfromtxt(p + '_vdb_faces.txt',dtype=np.int64)
        vdb_faces = vdb_faces - 1
        
        ivas_verts = np.genfromtxt(p + '_ivas_vertices.txt')
        ivas_faces = np.genfromtxt(p + '_ivas_faces.txt',dtype=np.int64)
        ivas_faces = ivas_faces - 1
        
        #direction cosines
        vdb_normals = np.genfromtxt('vertex_normals_vdb_' + filenames[i] + '.txt')
        alpha_rad = np.arccos(vdb_normals[:,0])
        beta_rad = np.arccos(vdb_normals[:,1])
        gamma_rad = np.arccos(vdb_normals[:,2])
        
        vdb_alpha = np.rad2deg(alpha_rad)
        vdb_beta = np.rad2deg(beta_rad)
        vdb_gamma = np.rad2deg(gamma_rad)
        
        ivas_normals = np.genfromtxt('vertex_normals_ivas_' + filenames[i] + '.txt')
        alpha_rad = np.arccos(ivas_normals[:,0])
        beta_rad = np.arccos(ivas_normals[:,1])
        gamma_rad = np.arccos(ivas_normals[:,2])
        
        ivas_alpha = np.rad2deg(alpha_rad)
        ivas_beta = np.rad2deg(beta_rad)
        ivas_gamma = np.rad2deg(gamma_rad)
        
        #Geometrically, the scalar triple product
        #
        #    a ⋅ ( b × c ) {\displaystyle \mathbf {a} \cdot (\mathbf {b} \times \mathbf {c} )} \mathbf{a}\cdot(\mathbf{b}\times \mathbf{c}) 
        #
        #is the (signed) volume of the parallelepiped def
        #ined by the three vectors given. Here, the parentheses may be omitted without causing ambiguity, since the dot product canno
        #t be evaluated first. If it were, it would leave the cross product of a scalar and a vector, which is not defined.
    
        # scalar triple product
        vdb_scalar_triple_products = []
        for i,f in enumerate(vdb_faces.tolist()):
        #        if i == 2:
        #            print(f)
        #            print(vdb_normals[f[0],:])
        #            print(vdb_normals[f[1],:])
        #            print(vdb_normals[f[2],:])
        #            print(np.cross(vdb_normals[f[1],:], vdb_normals[f[2],:]))
        #            bc_cross_product = np.cross(vdb_normals[f[1],:], vdb_normals[f[2],:])
        #            print(np.dot(bc_cross_product, vdb_normals[f[0],:]) )
                
                # a dot (b cross c)
        #        bc_cross_product = np.cross(vdb_normals[f[1],:], vdb_normals[f[2],:])
        #        scalar_triple_product = np.dot(bc_cross_product, vdb_normals[f[0],:]) 
        #        vdb_scalar_triple_products.append(scalar_triple_product)
                
                # or easier with the determinant with np.linalg.det(a)
                scalar_triple_product = np.linalg.det(np.array([vdb_normals[f[0],:],vdb_normals[f[1],:],vdb_normals[f[2],:]])) 
                vdb_scalar_triple_products.append(scalar_triple_product)
                
                
                
        ivas_scalar_triple_products = []
        for i,f in enumerate(ivas_faces.tolist()):
                
                # a dot (b cross c)
                bc_cross_product = np.cross(ivas_normals[f[1],:], ivas_normals[f[2],:])
                scalar_triple_product = np.dot(bc_cross_product, ivas_normals[f[0],:]) 
                ivas_scalar_triple_products.append(scalar_triple_product)
        
        vdb_color = 'gray'
        ivas_color = 'b'
                
        low_per = 10
        high_per = 90
        vdb_perc_low = np.percentile(np.abs(vdb_scalar_triple_products), low_per)
        vdb_perc_high = np.percentile(np.abs(vdb_scalar_triple_products), high_per)
        ivas_perc_low = np.percentile(np.abs(ivas_scalar_triple_products), low_per)
        ivas_perc_high = np.percentile(np.abs(ivas_scalar_triple_products), high_per)
        
        f, (ax2, ax1) = pl.subplots(1, 2, sharey=True,  sharex=True)
        
        if p == 'real':
            ##real
            pl.ylim(ymin=0,ymax=0.8)
            ax2.set_ylim(0,0.8)
            ax2.set_xlim(-20,450)
        elif p == 'dist':
        #dist
            ax2.set_ylim(-0.01,0.2)
            ax2.set_xlim(-20,650)
        elif p == 'sharp':
        ##sharp
            ax2.set_ylim(-0.01,0.14)
            ax2.set_xlim(-20,650)
            
        
        ax2.scatter(np.arange(len(ivas_scalar_triple_products)),np.abs(ivas_scalar_triple_products), marker=',', c=ivas_color , label= 'MC')        
        fontname = "Times New Roman"
        labels = ax2.get_xticklabels() + ax2.get_yticklabels()
        [label.set_fontname(fontname) for label in labels]
        

        
        
        ax2.axhline(y=ivas_perc_low,xmin=0,xmax=3,c=ivas_color,linewidth=2,zorder=0)
        
        
        ax2.axhline(y=ivas_perc_high,xmin=0,xmax=3,c=ivas_color,linewidth=2,zorder=0)
        
        ax2.tick_params(axis='both', which='major', labelsize=19)

        pl.xticks(fontname = "Times New Roman")  # This argument will change the font.
        pl.yticks(fontname = "Times New Roman")  # This argument will change the font.



        if p == 'real':
            ##real
            xticks = np.array([0,200,400])
            
        
        else:
            xticks = np.array([0,300,600])
            
        pl.xticks(xticks, fontname = "Times New Roman")  # This argument will change the font.
        pl.yticks(fontname = "Times New Roman")  # This argument will change the font.
        ax2.set_xlabel('triangle index ' ,fontsize=22, fontname = "Times New Roman")
        ax2.set_ylabel('|scalar triple product| / $nm^3$', fontsize=22, fontname = "Times New Roman")
        ax2.grid()

        
        pl.rcParams['legend.fontsize'] = 14.0
        ax2.legend(scatterpoints=1)
        
        
        
        ax1.scatter(np.arange(len(vdb_scalar_triple_products)),np.abs(vdb_scalar_triple_products),  marker='o', c=vdb_color, label= 'DMC')
        fontname = "Times New Roman"
        labels = ax1.get_xticklabels() + ax1.get_yticklabels()
        [label.set_fontname(fontname) for label in labels]        
        ax1.axhline(y=vdb_perc_low,xmin=0,xmax=3,c=vdb_color,linewidth=2,zorder=0)
    
        
       
        ax1.axhline(y=vdb_perc_high,xmin=0,xmax=3,c=vdb_color,linewidth=2,zorder=0)
        
        if p == 'real':
            ##real
            pl.ylim(ymin=0,ymax=0.8)
            ax1.set_ylim(0,0.8)
            ax1.set_xlim(-20,450)
        elif p == 'dist':
        #dist
            ax1.set_ylim(-0.01,0.2)
            ax1.set_xlim(-20,650)
        elif p == 'sharp':
        ##sharp
            ax1.set_ylim(-0.01,0.14)
            ax1.set_xlim(-20,650)
            
        if p == 'erbco':
            ##real
            xticks = np.array([0,200,400])
            
        
        else:
            xticks = np.array([0,300,600])
        
        ax1.tick_params(axis='x', which='major', labelsize=19)
        pl.yticks(fontname = "Times New Roman")  # This argument will change the font.
        pl.xticks(fontname = "Times New Roman")  # This argument will change the font.
        ax1.set_xlabel('triangle index ' ,fontsize=22, fontname = "Times New Roman")
#        ax1.set_ylabel('|scalar triple product| / $nm^3$', fontsize=22, fontname = "Times New Roman")
        ax1.grid()
        pl.rcParams['legend.fontsize'] = 14.0
        ax1.legend(scatterpoints=1)
        


            


        f.tight_layout()
        f.subplots_adjust(hspace=0)
    
        #f.savefig('/home/lukas/master_thesis_newdesign/figures/combined_ivas_vdb_' + p + '_stp_1000dpi.eps', format='eps', dpi=1000)
        
    
    #
    
    #fig2.savefig('/home/lukas/master_thesis_newdesign/figures/ivas_' + p + '_stp_1000dpi.eps', format='eps', dpi=1000)
