MATLAB CODE for reproducing the RESULTS OF FIG 5 of the paper, titled 

G. Tzimiropoulos, S. Zafeiriou and M. Pantic, "Robust and Efficient Parametric Face Alignment", ICCV 2011.

http://ibug.doc.ic.ac.uk/media/uploads/documents/iccv_final.pdf

@inproceedings{tzimiroREPFAICCV2011,
    author = {G. Tzimiropoulos and S. Zafeiriou and M. Pantic},
    pages = {1847--1854},
    booktitle = {Proceedings of IEEE Int’l Conf. on Computer Vision (ICCV 2011) },
    month = {November},
    title = {Robust and Efficient Parametric Face Alignment},
    year = {2011},
}

PLEASE CITE OUR PAPER ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Unless otherwise stated, all code written by G. Tzimiropoulos, S. Zafeiriou and M. Pantic 

Intelligent Behaviour Understanding Group (IBUG), Department of Computing, Imperial College London

Version: 1.0, 03/01/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yale database : http://cvc.yale.edu/projects/yalefacesB/yalefacesB.html

Code tested on 64-bit Windows 7 machine running MATLAB 7.11.0 (R2010b)

Please set the ICCV2011 folder as the current Matlab Folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN SCRIPTS:
main_Yale : reproduces the experiments (FIG 5) on the Yale Database.

The methods tested (Please see section 5.2 of the paper) are listed in the variable alg_list:
'affine_ic'             		corresponds to ref. [5] of the paper
'affine_ECC_ic'         		corresponds to ref. [14] of the paper 
'affine_ic_irls'        		corresponds to ref. [6] of the paper 
'affine_GaborFourier_ic'                corresponds to ref. [13] of the paper
'affine_GradientImages_ic'              corresponds to the method "GradientImages" also described in the paper (Please see section 5.2 of the paper)
'affine_GradientCorr_ic'                corresponds to the proposed method "GradientCorr"

All methods are tested in the framework of
Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
http://www.ri.cmu.edu/research_project_detail.html?project_id=515&menu_id=261

Set verbose = 1 (default) to show the fitting procedure. If verbose = 0, the fitting is not shown. 

For each script, the average frequency of convergence results are stored in the "results" variable. 
This is a 4-D matrix: results(subj, i, s, l)
subj denotes the subject under examination (10 subjects in total for Yale, 27 subjects in total for AR)
i    denotes the image pair under examination (10 image pairs for each subject for Yale, 3 image pairs for each subject for AR)
s    denotes the initial point RMS (parameter sigma in the paper, 100 tests are carried out for each s)
l    denotes the method examined (6 methods used) 
  l=1 corresponds to 'affine_ic'  
  l=2 corresponds to 'affine_ECC_ic'    
  l=3 corresponds to 'affine_ic_irls' 
  l=4 corresponds to 'affine_GaborFourier_ic' 
  l=5 corresponds to 'affine_GradientImages_ic'
  l=6 corresponds to the proposed 'affine_GradientCorr_ic'  
Example: results(3, 1, 10, 6) gives the average frequency of convergence (over 100 tests for initial point RMS = 10
of the proposed GradientCorr method tested on the 1st image pair of the 3rd subject.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FOLDERS:
1. Data    : the data used.
2. Methods : the implementations of the methods used. Please see above and section 5.2 of the paper.
             Note: the implementation of the proposed  "GradientCorr" is currently p-encoded. 
3. Tools   : common tools used for the implementation of Methods.
4. Results : scripts which plot the obtained results. (reproduce Fig 5 the paper)
             plot_results_Yale : reproduces Fig 5 a of the paper
             plot_results_AR   : reproduces Fig 5 b and c of the paper

 


