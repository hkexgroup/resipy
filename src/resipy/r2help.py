#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:46:04 2018

@author: jkl
"""

r2help = {
'flux_type': "<p>Define the flux type to be used (either 2D or 3D)</p>",
'res_matrix' : "<p><code>res_matrix</code> is 1 if a 'sensitivity matrix' is required for the converged solution (default). This matrix is not the Jacobian but is the diagonal of [J^T.W^T.W.J] which gives an idea of the mesh sensitivity (see equation 5.20 of Binley and Kemna, 2005)</p><p>Set <code>res_matrix</code> to 2 if the true resolution matrix is computed for a converged solution and the diagonal is stored (see equation 5.18 of Binley and Kemna,  2005), note that this calculation is more time consuming than the ‘sensitivity matrix’ option.</p><p>Set <code>res_matrix</code> to 3 if the sensitivity map is required and an output of the jacobian matrix and roughness matrix.</p><p>If neither sensitivity map or resolution matrix is required then set <code>res_matrix</code> to 0.</p>",
'singular_type' : "<p><b>Warning</b>: Note that singularity removal can only be applied is (a) the boundaries are infinite (far away from the region of interest) and (b) the z=0 plane defines the upper boundary (so a flat mesh). This is because the singularity remove uses the analytical solution to remove infinite potential around the electrode. If these two conditions applied, removing the singularity will make the forward computation more accurate.</p>",
'inverse_type' : "<p><code>inverse_type</code> is 0 for pseudo-Marquardt solution or 1 for regularised solution with linear filter (default) or 2 for regularised type with quadratic filter or 3 for qualitative solution or 4 for blocked linear regularised type.</p>",
'target_decrease': "<p><code>target_decrease</code> is a real number which allows the user to specify the relative reduction of misfit in each iteration. A value of 0.25 will mean that R2 will aim to drop the misfit by 25% (and no more) of the value at the start of the iteration. This allows a slower progression of the inversion, which can often result in a better convergence. If you set <code>target_decrease</code> to 0.0 then R2 will try to achieve the maximum reduction in misfit in the iteration (default).</p> ",
'qual_code' : "<p><code>qual_ratio</code> is 0 for qualitative comparison with forward solution, i.e. only when one observed data set is available, or <code>qual_ratio</code> is 1 if the observed data in <code>protocol.dat</code> contains a ratio of two datasets.</p>",
'scale' : "<p><code>scale</code> is a scaling factor for the mesh coordinates. This is usually 1.0 but if a standardised mesh is used, say for a unit circle, then this scaling factor is useful to adjust the mesh for a specific problem. Set <code>scale</code>=1 if you do not wish to change the coordinates of the mesh defined in <code>mesh.dat</code> (default)</p>",
'num_regions' : "<p><code>num_regions</code> is number of resistivity regions that will be specified either as starting condition for inverse solution or actual model for forward solution. The term “region” has no significance in the inversion – it is just a means of inputting a non-uniform resistivity as a starting model for inversion or for forward calculation.</p>",
'patch' : "<p><code>patch_size_x</code> and <code>patch_size_z</code> are the parameter block sizes in the x and z direction, respectively. We differentiate between parameter size and element size to allow faster computation. The larger the patch size the few parameters and the faster the inversion, however, if we increase it too much we will reduce the flexibility to create variation of resistivity. If computational time is not a problem then use a patch size of 1 for x and z (default). Note that the number of elements in the x direction must be perfectly divisible by <code>patch_size_x</code> and the number of elements in the z direction must be perfectly divisible by <code>patch_size_z</code> otherwise set them both to zero.</p>",
'data_type' : "<p><code>data_type</code> is 0 for true data based inversion or 1 for log data based (default to 1.0 for IP). Note that the latter should improve convergence but may not work for internal electrodes (e.g. borehole type) where the polarity can change due to resistivity distributions</p>",
'reg_mode' : "<p><code>reg_mode</code> is 0 for normal regularisation; or 1 if you want to include regularisation relative to your starting resistivity; or 2 if you wish to regularise relative to a previous dataset using the “Difference inversion” of LaBrecque and Yang (2000). If you select <code>reg_mode</code>=1 then you can define a regularisation parameter <code>alpha_s</code>. Note that if you select <code>reg_mode</code>=2 then <code>data_type</code> is automatically set to 0 irrespective of what is selected in the UI.</p>",
'tolerance' : "<p><code>tolerance</code> is the desired misfit (usually 1.0)</p>",
'max_iterations' : "<p><code>max_iterations</code> is the maximum number of iterations.</p>",
'error_mod' : "<p><code>error_mod</code> is 0 if you wish to preserve the data weights, 1 or 2 if you wish the inversion to update the weights as the inversion progresses based on how good a fit each data point makes. <code>error_mod</code>=2 is recommended – this is a routine based on Morelli and LaBrecque (1996). Note that no weights will be increased.</p>",
'alpha_aniso' : "<p>The smoothing factor used in the code is alpha. <code>alpha_aniso</code> is the anisotropy of the smoothing factor, set <code>alpha_aniso</code> > 1 for smoother horizontal models, <code>alpha_aniso</code> < 1 for smoother vertical models, or <code>alpha_aniso</code>=1 for normal (isotropic) regularisation.</p>",
#'alpha_s' : "<p><code>alpha_s</code> is the regularisation to the starting model (if you set <code>reg_mode</code> = 1 in Line 21). Set <code>alpha_s</code> to a high value (e.g. 10) to highly penalise any departure from this starting model. Note that <code>alpha_s</code> will stay fixed – if you set it too high then R2 may not converge. R2.out will report the value of alpha used to regularise smoothing within the image – the regularisation relative to a reference model is additional to this. The user may find setting <code>alpha_s</code> useful as a comparison of inversions from two runs with difference reference models allows an assessment of the depth of investigation following the approach of Oldenburg and Li (1999).</p>",
'errorParam' : "<h3>For resistivity Only</h3>\
                <p>It is advisable to estimate <code>a_wgt</code> and <code>b_wgt</code> from error checks in the field data (ideally from reciprocal measurements - not measures of repeatability). Typically for surface data <code>a_wgt</code> will be about 0.01 ohms and <code>b_wgt</code> will be about 0.02 (roughly equivalent to 2% error). Note that if you select <code>data_type</code>=1 then although the resistance data are transformed into log apparent conductivities the <code>a_wgt</code> and <code>b_wgt</code> parameters should still reflect the variance of the resistance.</p>\
                <h3>For IP</h3>\
                <p><code>min_error</code> is the minimum magnitude error (this is to ensure that very low errors are not assigned and is only used if <code>a_wgt</code> and <code>b_wgt</code> are both zero), <code>a_wgt</code> and <code>b_wgt</code> are error variance model parameters: <code>a_wgt</code> is the relative error of magnitudes; <code>b_wgt</code> is absolute error of phase values (in mrad). It is advisable to estimate <code>a_wgt</code> and <code>b_wgt</code> from error checks in the field data (ideally from reciprocal measurements - not measures of repeatability). Typically for surface data <code>a_wgt</code> will be about 0.02 (equivalent to 2% error), <code>b_wgt</code> will be typically 2 mrad for good data, but could be much higher.</p>\
                <h3>For both</h3>\
                <p>Note also that you can select to include individual errors (e.g. from fitted error models) for each measurement in the data input file protocol.dat – to do this <code>a_wgt</code> and <code>b_wgt</code> should be set to 0.0.</p>",
'rho_max': " <code>rho_min</code> and <code>rho_max</code> are the minimum and maximum observed apparent resistivity to be used for inversion (use large extremes if you want all data to be used). If data are ignored by R2 because of the apparent resistivity limits then these will be reported individually in R2.log. NOTE: that the apparent resistivity calculations assume that you have set the ground surface to Y=0 and that the ground surface is flat.<p>",
'num_xy_poly' : "<p>where <code>num_xy_poly</code> is the number of x,y co-ordinates that define a polyline bounding the output volume. If <code>num_xy_poly</code> is set to zero then no bounding is done in the x-y plane. The co-ordinates of the bounding polyline follow in the next line. NOTE: the first and last pair of co-ordinates must be identical (to complete the polyline). So, for example, if you define a bounding square in x,y then you must have 5 co-ordinates on the polyline. The polyline must be defined as a series of co-ordinates in sequence, although the order can be clockwise or anti-clockwise (see examples later). NOTE: R2 stores the vertical co-ordinates for nodes in a quadrilateral mesh with a convention positive upwards. For example, if the ground surface has an elevation of 0m and you wish to output to a depth of 8m then y=-8m must be used for the lower boundary of the polygon. Similarly, if the ground surface elevation is 100m and you wish to output to a depth of 8m then y=-92m must be used for the lower boundary of the polygon. This was not the convention for v2.7a and so any input files created for that version must be changed (this only applies to line 26). If a triangular mesh is used then the co-ordinates specified in the mesh file are used and the above comments about sign convention do not apply. </p>",
'modErr' : "<p>Compute modelling error due to the mesh by doing a forward modelling with a homogeneous resistivity of 100 Ohm.m. This error will be geometrically added to the reciprocal error if the latest is available.</p>",
'parallel' : "<p>Run inversions in parallel according to the number of cores "
              "of your computer. Note that you won't be able to see the "
              "progress of the inversion. Don't "
              "do this for 3D surveys as the required memory might be huge.</p>",
'no_improve' : "<p> <code>no_improve</code> is termination criteria such that"
                " if during two iterations the misfit doesn\'t change by "
                "<code>no_improve</code> (%) then the inverse solution is stopped.</p>",
'inv_type3D' : "<p>0 for normal regularisation; 1 for background regularisation; 2 for difference regularisation. Note that option 2 is not a difference inversion (as per LaBrecque and Yang, 2001) but can be used for this purpose"
                "with a modification to the data in <code>protocol.dat</code>.</p>",
'alpha_s' : "<p><code>alpha_s</code> is an additional penalty factor applied to the starting resistivity"
            "model. If <code>alpha_s</code> is 1.0 then the regularisation applies the same"
            "weight to smoothing the model as to constraining to the background model. A smaller"
            "(no zero) value of <code>alpha_s</code> will retain some constraint to the background model.</p>",
'cginv_tolerance' : "<p><code>cginv_tolerance</code> is the tolerance (typically 0.0001) for the conjugate gradient "
                    "solution of the inverse equations.</p>",
'cginv_maxits' : "<p><code>cginv_maxits</code> is the maximum number of iterations for"
                "the conjugate gradient solution of the inverse problem (this could be set to a high value < "
                "number of parameters or could be set low, say 50, to gain an approximate solution). A "
                "value of 500 for <code>cginv_maxits</code> will achieve a satisfactory solution for most problems, "
                "however, setting the value to 50 will lead to a faster execution.</p>",
'alpha_max' : "<p>The regularisation (or smoothing) parameter, alpha, is optimised each iteration by carrying "
                "out a line search. <code>alpha_max</code> is maximum starting value of reqularised scalar and "
                "<code>num_alpha_steps</code> (usually 10) is the number of alpha values used each iteration to search "
                "for the optimum alpha. Set alpha_max to a large number (say 10e10) if you do not wish to "
                "limit the maximum value of alpha. If you wish to specify a minimum value for the starting "
                "alpha then set <code>alpha_max</code> to a negative value of the minimum starting value. For example, "
                "setting <code>alpha_max</code> to -20 will mean that the starting value of alpha is at least 20. An "
                "alternative solution approach is to use one alpha at each of the inverse iterations, i.e. there "
                "is no line search for an optimum at each iteration. This can help result in a smoother final "
                "model. If <code>num_alpha_steps</code> is set to 1 then <code>alpha_max</code> is the starting value of alpha. This "
                "value is then used in subsequent iterations unless the total objective function does not drop "
                "by 5% in a subsequent iteration, or if the data misfit increases. If this happens then alpha "
                "is reduced by 50% in the following iteration. If this approach is taken then it is advisable "
                "to set <code>max_iterations</code> to at least 20.</p>",
'min_step' : "<p><code>min_step</code> is the minimum step length for attempting to improve solution. This is "
            "usually set to 0.001 to 0.01.</p>",
'notCropping' : '<p>If checked, this will prevent the mesh to be cropped to the region of interest after inversion.',
'cropBelowFmd' : '<p>If checked, this will crop out the mesh (inverted plot) below fine/coarse boundary depth (set in mesh tab).<br><i>Only works on 2D meshes</i>.',
'modelDOI' : "<p>If checked, two more background constrained inversion will be run with the normal background and a background ten times more resistive."
            " The ratio of the difference of inverted value divided by the difference in initial resistivity gives"
            " an estimate of the depth of investigation of the inversion (Oldenburg and Li, 1999)</p>",
'sensDOI' : "<p>Sensitivity based DOI. If a dashed line will be drawn on the inverted section at 0.001*maxSensitivity."
            "Note that this is not the DOI based on the Oldenburg and Li method.</p>",
'txSign' : "<p>Check if the polarity of provided transfer resistanses is correct in each survey and if not, automatically corrects them."
            "<br><i>Note: this assumes a provided survey is on a flat 2D profile.</i></p>"
}

