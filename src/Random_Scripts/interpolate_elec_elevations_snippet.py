#interpolating electrode elevations snippet

k=R2(typ='R3t',dirname=baseline_dir)
k.setTitle('Hollin Hill reference')
#find closest peg survey to reference date, there is some sorting of survey dates 
days = timeDeltas(ref_date,peg_date)
date_idx = np.argmin(np.abs(days))
elec_x,elec_y,elec_z = gps.get_peg_pos(elec_df,date_idx) # start electrode positions 
elec = np.zeros((len(elec_x),3))
elec[:,0] = elec_x ; elec[:,1] = elec_y ; elec[:,2] = elec_z


topo_x,topo_y,topo_z = gps.get_peg_pos(surf_df,date_idx) # start surface positions
topo_x = np.array(topo_x);topo_y = np.array(topo_y);topo_z = np.array(topo_z) 


topo = np.zeros((len(topo_x[ikeep]),3))
topo[:,0] = topo_x 
topo[:,1] = topo_y 
topo[:,2] = topo_z
#for now interpolate electrode elevation due to 2m descrepancy 
elec[:,2] = interp2d(elec_x,elec_y,topo[:,0],topo[:,1],topo[:,2])