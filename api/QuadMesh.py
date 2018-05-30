#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 20:04:59 2018

@author: jkl
"""

import numpy as np


class QuadMesh(object):
    def __init__(elecPos=[], elec=[], spacing=1, nnode=4, patchx=1, patchy=1):
    ''' generate a quadrilateral mesh
        
    elec = electrode position as a numpy array (3 columns)
    nnode = number of nodes between two electrode, that control the
        mesh resolution
    elecPos = electrode number
    patchx, patchy : how groupd the element
    '''
    
    if len(elecPos) == 0:
        x = np.arange(len(elec), spacing)
        y = np.zeros(x)
        z = np.zeros(x)
        elecPos = np.vstack([x,y,z])
    
    nelec = len(elecPos)
    spacing = np.sqrt(np.diff(elecPos[:,2], axis=0)**2) # spacing between two first electrodes
    dx = spacing/nnode # regular dx, can implement an irregular dx
    dy = spacing/4

    itopo = np.sum(elecPos[:,1]) != 0 # TRUE if topo is specified
    meshx = np.zeros(nelec*nnode+4*nnode)
    meshx[2*nnode:-2*nnode] = np.arange(-nnode*dx, nelec*dx+nnode*dx, dx)
    meshx[:2*nnode] = 
    

''' javascript code for r2.html below
    
    
createRectangularMesh: function() {
    this.num_electrodes = this.getParameter("num_electrodes")
    this.elecSpacing = this.getParameter("elecSpacing")
    this.surveyLength = (this.num_electrodes-1)*this.elecSpacing
    this.patchx = this.getParameter("patchx")
    this.patchy = this.getParameter("patchy")
    var dx = this.elecSpacing/this.elemx // CHANGE HERE
    var dy = this.elecSpacing/4
    
    // for x
    var meshx = [0.0]
    var xnnodes = ((this.num_electrodes-1)*this.elecSpacing)/dx
    for (i=0; i<xnnodes+8; i++) {
      meshx.push(meshx[i]+dx)
    }
    console.log(meshx.length)
    for (i=0; i<8; i++) {
      meshx.unshift(meshx[0]-dx)
    }
    var dist = 8*dx+2*dx
    var negext2 = []
    var posext2 = []
    while (dist < 4.*this.surveyLength) {
      negext2.unshift(-dist)
      posext2.push(this.surveyLength+dist)
      dist = dist*2
    }
    // check with patchx and patchy
    var tmpLength = negext2.length+meshx.length+posext2.length
    var remain = (tmpLength-1) % this.patchx // -1 !!
    //console.log("initial remain="+remain)
    //console.log("initial length = " + tmpLength)
    if (remain == 0) {
      console.log("Number of X node perfectly match patchx.")
    }
    else if (remain == 1) {
      meshx = meshx.slice(1,meshx.length)
    }
    else if (remain == 2) {
      meshx = meshx.slice(1,meshx.length-1)
    }
    else if (remain == 3) {
      meshx = meshx.slice(2,meshx.length-1)
    }
    var tmpLength = negext2.length+meshx.length+posext2.length
    var remain = tmpLength % this.patchx

    meshx = negext2.concat(meshx.concat(posext2))
    //console.log(meshx)
    //console.log("latest element" + meshx[meshx.length-1])
    console.log("Total number of x node = " + meshx.length)

    // for y
    var meshy = [0.0]
    const f1 = 1.1
    const f2 = 1.5
    var dist = 0.05
    for (i = 0; i < 20; i++) {
      dist = dist*f1
      console.log(dist)
      meshy.push(meshy[meshy.length-1] + dist)
    }
    for (i = 0; i < 8; i++) {
      dist = dist*f2
      meshy.push(meshy[meshy.length-1] + dist)
    }
    console.log(meshy)
    console.log(meshy.length)

    var remain = meshy.length % this.patchy
    if (remain == 0) {
      console.log("Patch y perfectly match meshy lenth.")
    }
    else if (remain == 1) {
      meshy = meshy.slice(0,meshy.length-1)
    }
    else if (remain == 2) {
      meshy = meshy.slice(0,meshy.length-2)
    }
    else if (remain == 3) {
      meshy = meshy.slice(0,meshy.length-3)
      // warning :not so good this way -> add more elements ?
    }
    this.meshx = meshx
    this.meshy = meshy
    // TODO if patchx or patchy is bigger than 4 ?!

    this.injectRegions(document.getElementById("num_regions"))
    document.getElementById("meshOutput").innerHTML = "Mesh generated (numnp_x = "
    + this.meshx.length + " and numnp_y = " + this.meshy.length + ")."

    // create topo
    this.topo = []
    if (document.getElementById("topoCheck").checked) {
      // read the topo points specify
      x = document.getElementById("topo").value
      var tmp = x.split("\n") // read pair position topography
      var tmpx = []
      var tmpy = []
      for (i = 0; i < tmp.length; i++) {
        if (tmp[i].match(/\d+/)) {
          tmpx.push(parseFloat(tmp[i].split(/[ \t]+/)[0]))
          tmpy.push(parseFloat(tmp[i].split(/[ \t]+/)[1]))
        }
      }
      console.log("tmpx = " + tmpx)
      console.log("tmpy = " + tmpy)
      var tmpxmin = Math.min.apply(null, tmpx)
      var tmpxmax = Math.max.apply(null, tmpx)
      console.log("tmpminx = " + tmpxmin + " tmpxmax = " + tmpxmax)
      var newtopo = []
      for (i = 0; i < this.meshx.length; i++) {
        if (this.meshx[i] <= tmpxmin) {
          this.topo.push(tmpy[0])
          console.log("keep left")
        } else if (this.meshx[i] >= tmpxmax) {
          this.topo.push(tmpy[tmpy.length-1])
          console.log("keep right")
        } else {
          // interpolation
          var fstop = 0
          for (j = 0; j < tmp.length; j++) {
            //console.log("loop : j = " + j + " this.mesh[i] = " + this.meshx[i] + "; tmpx[j] = " + tmpx[j])
            if (this.meshx[i] < tmpx[j]) {
              fstop = 1
              break
            }
            if (this.meshx[i] == tmpx[j]) {
              fstop = 2
              break
            }
          }
          if (fstop == 1) {
            //console.log("j = " + j)
            var value = tmpy[j-1] + (tmpy[j]-tmpy[j-1]) / (tmpx[j]-tmpx[j-1]) * (this.meshx[i]-tmpx[j-1])
            //console.log("computed topo : " + value)
            this.topo.push(value)
          }
          else if (fstop == 2) {
            this.topo.push(tmpy[j])
          }
        }
      }
    } else {
      for (i = 0; i < this.meshx.length; i++) {
        this.topo.push(0.0)
      }
    }
    console.log("Generated topography : " + this.topo)

    // plot the mesh with plotly
    this.plotMesh()
  },
    
    
'''
    