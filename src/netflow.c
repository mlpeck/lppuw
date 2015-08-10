#include "lppuw.h"

/* network flow model for phase unwrapping */

void netflow(const int *ndx, 
           const int *ndy,
           int *nconstr,
           const int *k2x1, 
           const int *k2x2,
           const int *k2y1,
           const int *charge,
           double *objective,
           double *objval,
           double *soln,
           int *retval)
{
  lprec *lp;
  int nvars = 2*((*ndx) + (*ndy));
  int i, ret, colno[8];
  int constrtype = EQ;
  double row[8] = {-1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0};
  lp = make_lp(*nconstr, nvars);
  if (lp == NULL) {
  	*retval = MAKELP_FAIL;
  	return; /* couldn't construct a new model... */
  }
  set_verbose (lp, CRITICAL);
  set_minim(lp);
  set_add_rowmode(lp, true);

  if (objective == NULL) {
    objective = malloc((nvars+1) * sizeof(double));
    objective[0] = 0.0;
    i=1;
    while (i <= nvars) objective[i++] = 1.0;
  }
  ret = set_obj_fn(lp, objective);
  if (!ret) {
  	*retval = SETOBJ_FAIL;
  	return;
  }
  
  for (i = 0; i < (*nconstr); i++) {
  	colno[2] = (*ndx) + (colno[0] = k2x1[i]);
  	colno[3] = (*ndx) + (colno[1] = k2x2[i]);
  	colno[6] = (*ndy) + (colno[4] = k2y1[i] + 2*(*ndx));
  	colno[7] = (*ndy) + (colno[5] = colno[4] + 1);
/*
  	ret = add_constraintex(lp, 8, row, colno, constrtype, (double) charge[i]);
*/
	ret = set_rowex(lp, i+1, 8, row, colno);
  	if (!ret) {
  		*retval = ADDCONSTR_FAIL;
  		return;
  	}
  	set_rh(lp, i+1, (double) charge[i]);
  	set_constr_type(lp, i+1, constrtype);
  }
  set_add_rowmode(lp, false);
/* 
  set_basiscrash(lp, CRASH_MOSTFEASIBLE);
*/
  set_scaling(lp, SCALE_NONE);
  /* try adding presolve */
  /*
  set_presolve(lp, PRESOLVE_ROWS, get_presolveloops(lp));
  */
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

