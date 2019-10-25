##' Decision Rules
##'
##' Different rules for making decisions in the framework of belief functions
##' @export
##' @param mass The matrix containing the masses. Each column represents a piece of mass.
##' @param criterion The decision baseline:
##'
##'     criterion=1 maximum of the plausibility
##'
##'		criterion=2 maximum of the credibility
##'
##'		criterion=3 maximum of the credibility with rejection
##'
##'		criterion=4 maximum of the pignistic probability
##'
##'		criterion=5 Appriou criterion (decision onto \eqn{2^\Theta})
##'
##'		criterion=6 Distance criterion (decision onto a given subset (sDec) of \eqn{2^\Theta})
##'
##' @param r The parameter in BayesianMass function. If criterion 5 is used, it should be given. 
##' Otherwise it will be set to the default value 0.5.
##' @param sDec The parameter for the set on which we want to decide. It should be a subset of {1,2,...,\eqn{2^n}}, where \eqn{n} is the number of elements in \eqn{\Theta}. If criterion 6 is used, it should be given; Otherwise it will be set as the default value \eqn{2^\Theta} 
##' @param D The parameter for the used matrix in Jousselme distance. If criterion 6 is used, it should be given. Otherwise it will be  set as default
##' Otherwise it will be calculated.
##' @return The decision vector. E.g., in classification problem, class labels. 
##' @examples
##' m1=c(0,0.4, 0.1, 0.2, 0.2, 0, 0, 0.1);
##' m2=c(0,0.2, 0.3, 0.1, 0.1, 0, 0.2, 0.1);
##' m3=c(0.1,0.2, 0, 0.1, 0.1, 0.1, 0, 0.3);
##' 
##' m3d=discounting(m3,0.95);
##' 
##' M_comb_Smets=DST(cbind(m1,m2,m3d),1);
##' M_comb_PCR6=DST(cbind(m1,m2),8);
##' 
##' class_fusion=decisionDST(M_comb_Smets,1)
##' class_fusion=decisionDST(M_comb_PCR6,1)
##' class_fusion=decisionDST(M_comb_Smets,5,0.5)
##' class_fusion=decisionDST(cbind(M_comb_Smets,M_comb_PCR6),1)
##' sDec<-c(2,3,4)
##' class_fusion=decisionDST(M_comb_Smets,6, sDec = sDec)
##' 
decisionDST <- function (mass, criterion, r = 0.5, sDec = 1:nrow(mass), D = Dcalculus(nrow(mass))){

  if (is.vector(mass) || (is.matrix(mass) && nrow(mass) == 1)) {
    mass = matrix(mass,, 1)
  }
	lm=nrow(mass);
	nbvec_test=ncol(mass);
	nbclasses=round(log2(lm));



	class_fusion=c();
	for(k in 1:nbvec_test){
		masstmp=mass[,k];

		if(criterion==1){
			# case 1
			plau=mtopl(masstmp);
			ii=1:nbclasses;
			plau_singl=plau[1+2^(ii-1)];
			indice=which.max(plau_singl);
			class_fusion=c(class_fusion,indice);
		}else if(criterion==2||criterion==3){
			# case {2, 3}
      # browser()
			croy=mtobel(masstmp);
			ii=1:nbclasses;
			croy_singl=croy[1+2^(ii-1)];
	        valuemax=max(croy_singl);
			indice=which.max(croy_singl);
				if(criterion==3){
					indice_complementaire=0;
					for (i in seq(nbclasses,indice,by=-1)){
						indice_complementaire=indice_complementaire+2^(nbclasses-(nbclasses-i+1));
				    }
			    if (valuemax>=croy[indice_complementaire]){
						class_fusion=c(class_fusion,indice);
				}else{
						class_fusion=c(class_fusion,0);
				}
				}else{
				  class_fusion=c(class_fusion,indice);
        }
		}else if(criterion==4){
			# case 4
				pign=mtobetp(t(masstmp));
				indice=which.max(pign);
				class_fusion=c(class_fusion,indice);
		}else if(criterion==5){
			# case 5
				plau=mtopl(t(masstmp));
				lambda=1;
				md=BayesianMass(lambda,r,nbclasses);
				indice=which.max(plau*md);
				class_fusion=c(class_fusion,indice);
		}else if(criterion==6){
			# case 6
			   sizeSD <- length(sDec)
			   distJ <- c()
			   for (i in 1:sizeSD){
			 	 mSD <- matrix(0,lm,1)
				 mSD[sDec[i]] <- 1
				 distJ <- c(distJ, JousselmeDistance(masstmp, mSD, D))
		       }
		       indice=which.min(distJ)
		       class_fusion=c(class_fusion,sDec[indice])
		 }else{
				stop('ACCIDENT: The critertion given is not right\n')
		}
  }
	return(class_fusion)
}

Dcalculus <- function(lm) {
  # computing the table of conflict for natoms = round(log2(lm)) classes
  
  natoms = round(log2(lm))
  ind = list() 
  if (2^natoms == lm) {
    ind[[1]] = c(1)
    ind[[2]] = c(2)
    step = 3
    while (step < 2^natoms) {
      ind[step] = step
      step = step + 1
      indatom = step
      for (step2 in 2: (indatom - 2)) {
        ind[[step]] = sort(union(ind[[indatom - 1]], ind[[step2]]));
        step = step + 1
      }
    } 
    out = matrix(0, 2^natoms, 2^natoms)
    
    for (i in 1:2^natoms) {
      for (j in 1:2^natoms) {
        out[i, j] = length(intersect(ind[[i]], ind[[j]]))/length(union(ind[[i]], ind[[j]]))
      }
    }
  } else{
    stop("ACCIDENT in Dcalculus: length of input vector not OK: should be a power of 2")
  }
  return(out)
}

JousselmeDistance <- function(m1, m2, Tjaccard) {
  if (length(m1) == length(m2)) {
    if(missing(Tjaccard)){
      Tjaccard = Dcalculus(length(m1))
    }
    m_diff = matrix(m1 - m2, length(m1) ,1)
    out = sqrt(t(m_diff) %*% Tjaccard  %*% m_diff/2)
  } else {
    stop("ACCIDENT in JousselmeDistance: the size of the both masses m1 and m2 is different\n")
  }
  return(out)
} 




