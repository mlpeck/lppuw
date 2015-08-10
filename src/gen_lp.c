
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "lp_lib.h"

/* generic lp definition */

void gen_lp(const int *dir, 
           int *nvars, 
           double *objective, 
           int *nconstr,
           double *constr,
           int *constrtype,
           double *rhs,
           double *objval,
           double *soln)
{
  lprec *lp;
  int i, j, ret, *colno;
  double *row = NULL;

  lp = make_lp(0, *nvars);
  if (lp == NULL)
    return; /* couldn't construct a new model... */

  set_verbose (lp, CRITICAL);
  ret = set_obj_fn(lp, objective);
  if (ret == 0) return;
  
  if (*dir == 1)
    set_maxim(lp);
  else
    set_minim(lp);
  
  /* create space large enough for one row */
  colno = (int *) malloc((*nvars) * sizeof(*colno));
  row = (double *) malloc((*nvars) * sizeof(*row));
  if((colno == NULL) || (row == NULL))
    return;

  set_add_rowmode(lp, true);  /* makes building the model faster if it is done rows by row */
  
  for (i = 0; i < *nvars; i++) colno[i] = i+1;
  
  for (i = 0; i < *nconstr; i++) {
  	for (j = 0; j < *nvars; j++)
  		row[j] = *(constr + i + (*nconstr)*j);
  	ret = add_constraintex(lp, *nvars, row, colno, constrtype[i], rhs[i]);
  }

  set_add_rowmode(lp, false); /* rowmode should be turned off again when done building the model */
  print_lp(lp);

  ret = solve(lp);
  if (ret != 0) {
  	delete_lp(lp);
  	return;
  }
  *objval = get_objective(lp);
  get_variables(lp, soln);
  delete_lp(lp);
  return;

}

