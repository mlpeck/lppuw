#include "lppuw.h"
#include <stdio.h>

extern int ind(int i, int j, int k);

/* rows, cols, contain the indexes of each block, ordered by block number */

void jigsudoku(const int *givens, const int *ngiven, 
	int *rows, int *cols, double *soln, int *retval) {
	lprec *lp;
	int ret, i = 1, j, k, l, block, ib, nc = 1, gc, colno[9];
	int nconstr = 324+9*(*ngiven), nvars=729;
	double rhs = 1.0, row[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
	
	lp = make_lp(nconstr, nvars);
	set_verbose(lp, CRITICAL);
	ret = set_obj_fnex(lp, 1, &rhs, &i);
	set_add_rowmode(lp, true);
	l = 0;
	/* set the given values */
	
	for (j=0; j<9; j++) {
		for (i=0; i<9; i++) {
			if ((k = givens[l++]) > 0) {
				for (ib=0; ib<9; ib++) {
				  gc = ind(i, j, ib);
				  set_rowex(lp, nc, 1, row, &gc);
				  if (ib == k-1) rhs=1.0; else rhs = 0.0;
				  set_rh(lp, nc, rhs);
				  set_constr_type(lp, nc++, EQ);
				}
 			}
		}
	}
	
	/* Each cell has one entry */
	
	for (j=0; j<9; j++) {
		for (i=0; i<9; i++) {
			for (k=0; k<9; k++) colno[k] = ind(i, j, k);
			set_rowex(lp, nc++, 9, row, colno);
		}
	}
	
	/* Each row gets distinct numbers */
	
	for (i=0; i<9; i++) {
		for (k=0; k<9; k++) {
			for (j=0; j<9; j++) colno[j] = ind(i,j,k);
			set_rowex(lp, nc++, 9, row, colno);
  		}
	}
	
	/* Each column gets distinct numbers */
	
	for (j=0; j<9; j++) {
		for (k=0; k<9; k++) {
			for (i=0; i<9; i++) colno[i] = ind(i,j,k);
			set_rowex(lp, nc++, 9, row, colno);
		}
	}
	
	/* Each block gets distinct numbers */
	
	for (k=0; k<9; k++) {
	  l = 0;
	  for (block=0; block < 9; block++) {
	    for (i = 0; i < 9; i++) {
	      colno[i] = ind(rows[l], cols[l], k);
	      l++;
	    }
	    set_rowex(lp, nc++, 9, row, colno);
	  }
	}
	
	/* set constraint and rhs values for all rows */
	
	rhs = 1.0;
	for (i=1+9*(*ngiven); i<nc; i++) {
		set_rh(lp, i, rhs);
		set_constr_type(lp, i, EQ);
	}
	set_add_rowmode(lp, false);
	for (i=1; i<=nvars; i++) set_binary(lp, i, true);
	ret = solve(lp);
	if (ret != 0) {
  		delete_lp(lp);
  		*retval = ret;
  		return;
  	}
  	get_variables(lp, soln);

  	delete_lp(lp);
  	*retval = 0;
	return;
}
