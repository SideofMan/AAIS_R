# AAIS_R
R code for Bin Lau's Adaptive Annealed Importance Sampling algorithm

Object descriptions:

Proposal: has 4 parts  
	M (number of components)  
	W (weights vector of the components)  
	Mu (mean vectors of the components)  
	Sigma (covariance matrices of the components)  
X: has 13 parts  
	Values (points)  
	N (number of points)  
	Resp (output of t-distribution of points for all components)  
	Proposal (combination of Resp's with component weights)  
	logProposal (log of Proposal)  
	logAnnealTarget (eqn 7 in paper)  
	logWeight (log of weights of each point)  
	NormalizedWeight (normalized weight of each point)  
	Root (the component it was sampled from)  
	logProposal0 (log of original Proposal)  
	logPrior (log of Prior points)  
	logLike (log of likelihood)  
	logTarget (log of Target function, log of product between Prior and likelihood)  

Function descriptions:  

main_AAIS:  
	The main script that runs the whole algorithm.  
	It takes in a target function and runs the AAIS algorithm, outputting a Student's t-mixture distribution that approximates the target function  

Merge:  
	Inputs: Proposal object, X object, Cor_thr, df, dim  
	Outputs: Updated Proposal objectm and updated X object  
	Merges components (section 3.2 in the paper)  

Lightspeed helper functions:  
	These functions help speed up parts of the algorithm  

ISEM_astro:  
	Inputs: X object, Proposal object, Proposal0 object (initial Proposal object), df, Target_function, gama (annealing schedule), j (gama iteration number), data  
	Outputs: Updated Proposal objectm and updated X object  
	Performs one iteration of EM with importance sampling

log_t_pdf:
	Inputs: X points, Mu (mean vector), Sigma (covariance matrix), df
	Outputs: log of student-t's distribution of X points
	The log of a student's t-distribution function

t_mix_sample:  
	Inputs: Proposal object, X object, df  
	Outputs: X object with N sampled points  
	Samples N points from the Student's t-mixture distribution, it calls function t_mixture_sampling  

t_mixture_sampling:  
	Inputs: M (number of components), W_m (weights of components), Mu (mean vectors), Sigma (covariance matrices)  
	Outputs: one sampled point and the index of the component it was sampled from  
	Samples 1 point from the Student's t-mixture distribution  

t_mix_update_v2:  
	Inputs: X object, Proposal object, df  
	Outputs: Updated Proposal object  
	Initial EM algorithm to adjust proposal model  

t_mixture_pdf:  
	Inputs: X points, M (number of components), W (weights of components), Mu (mean vectors), Sigma (covariance matrices), df  
	Outputs: t-mixture distribution output of X points  
	pdf of the Student's t-mixture distribution  

t_pdf:  
	student's t-distribution function  
