#ifndef ERGMM_STRUCTS_H
#define ERGMM_STRUCTS_H

#define PROP_NONE (65535-1)
#define PROP_ALL (65535-2)

/* The structure to house the state of MCMC for a given iteration.
 * Includes only values that change from iteration to iteration.
 * Also includes the "sub-states" of MCMC: proposed Z and proposed coef,
 * as well as cached values for latent distances.
 * The variables *_llk_old are to store the respective variables from the
 * previous call to the likelihood. The variables *_old (no "llk") are to store
 * the original value of MH-updated variables when they are being proposed.
 */

typedef struct {
  double **Z, *coef, **Z_mean, *Z_var, *Z_pK;
  unsigned int *Z_K;
  double llk, **lpedge, lpZ, lpLV, lpcoef;
} ERGMM_MCMC_Par;

typedef struct {
  ERGMM_MCMC_Par *state,*prop;
  double **Z_bar,*deltas, *pK;
  unsigned int *n;
  unsigned int prop_Z, prop_coef, prop_LV, after_Gibbs;
  unsigned int *update_order;
} ERGMM_MCMC_MCMCState;

/* The structure to house the settings of MCMC: constants that, while
 * they affect the sampling, are not a part of the posterior distribution.
 */
typedef struct {
  double Z_delta;
  double **group_deltas;
  unsigned int group_prop_size;
  unsigned int sample_size, interval;
  unsigned int accept_all; // debugging option: accept all MH proposals
} ERGMM_MCMC_MCMCSettings;

/* The structure to house the parameters of the prior distribution. 
 */
typedef struct {
  double Z_mean_var, Z_var, Z_var_df, *coef_mean, *coef_var, Z_pK;
} ERGMM_MCMC_Priors;

/* The structure to house the MCMC draws.
 */
typedef struct {
  double *llk, *lpZ, *lpcoef, *lpLV;
  double *Z, *Z_rate_move, *coef, *coef_rate, *Z_mean, *Z_var, *Z_pK;
  int *Z_K;
} ERGMM_MCMC_ROutput;


/* The structure to house the data on which the posterior is
   conditioned, and model choice.
*/
struct ERGMM_MCMC_Model_struct{
  unsigned int dir;
  int **iY;
  double **dY;
  double ***X;
  unsigned int **observed_ties;

  double (*lp_edge)(struct ERGMM_MCMC_Model_struct*,ERGMM_MCMC_Par*,unsigned int,unsigned int);
  double (*E_edge)(struct ERGMM_MCMC_Model_struct*,ERGMM_MCMC_Par*,unsigned int,unsigned int);

  double lp_Yconst;
  int *iconst;
  double *dconst;
  unsigned int verts, latent, coef, clusters;
} ;

typedef struct ERGMM_MCMC_Model_struct ERGMM_MCMC_Model;

#endif
