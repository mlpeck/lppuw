#include "lppuw.h"

#define DMAX(a, b) ((a)>=(b) ? (a):(b))

/* modified assignment problem as used by brcutpuw */

void mod_assign(const int *ncp, 
           const int *ncm,
           double *objective, 
           double *objval,
           double *soln,
           int *retval)
{
  lprec *lp;
  int nvars = (*ncp)*(*ncm) + (*ncp) + (*ncm);
  int nconstr = (*ncp) + (*ncm);
  int i, j, nr=1, ret, maxel = DMAX(*ncp, *ncm) + 1, *colno;
  int constrtype = EQ;
  double rhs = 1.0;
  double *row = NULL;

  lp = make_lp(0, nvars);
  if (lp == NULL) {
  	*retval = MAKELP_FAIL;
  	return; /* couldn't construct a new model... */
  }

  set_verbose (lp, CRITICAL);
  set_minim(lp);
  
  ret = set_obj_fn(lp, objective);
  if (!ret) {
  	*retval = SETOBJ_FAIL;
  	return;
  }
  set_add_rowmode(lp, true);
  
  /* create space large enough for one row */
  colno = (int *) malloc(maxel * sizeof(int));
  row = (double *) malloc(maxel * sizeof(double));
  if((colno == NULL) || (row == NULL)) {
    *retval = ADDCONSTR_FAIL;
    return;
  }
  
  for (i = 0; i < maxel ; i++) row[i] = 1.0;
  
  for (i = 0; i < *ncp; i++) {
  	for (j = 0; j < *ncm; j++)
  		colno[j] = i + j*(*ncp) + 1;
  	colno[j] = i + (*ncp)*(*ncm) + 1;
  	ret = add_constraintex(lp, *ncm+1, row, colno, constrtype, rhs);
/*	ret = set_rowex(lp, nr, *ncm+1, row, colno); */
  	if (!ret) {
  		*retval = ADDCONSTR_FAIL;
  		return;
  	}
/*
  	set_rh(lp, nr, rhs);
  	set_constr_type(lp, nr++, constrtype);
*/
  }
  
  for (i = 0; i < *ncm; i++) {
  	for (j = 0; j < *ncp; j++)
  		 colno[j] = i * (*ncp) + j + 1;
  	colno[j] = (*ncp)*(*ncm) + (*ncp) + i + 1;
  	ret = add_constraintex(lp, *ncp+1, row, colno, constrtype, rhs);
/*	ret = set_rowex(lp, nr, *ncp+1, row, colno); */
  	if (!ret) {
  		*retval = ADDCONSTR_FAIL;
  		return;
  	}
/*
   	set_rh(lp, nr, rhs);
  	set_constr_type(lp, nr++, constrtype);
*/
 }

  set_add_rowmode(lp, false);
  
/*  set_basiscrash(lp, CRASH_MOSTFEASIBLE); */

  ret = solve(lp);
  if (ret != 0) {
  	delete_lp(lp);
  	*retval = ret;
  	return;
  }
  *objval = get_objective(lp);
  get_variables(lp, soln);
  delete_lp(lp);
  *retval = 0;
  return;

}

