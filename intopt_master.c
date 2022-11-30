#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <C:/Users/roosb/intopt/Lab4/labs_header.h>

double epsilon = 1e-6; 
typedef struct node_t node_t;
typedef struct simplex_t simplex_t;
typedef struct Set Set;

struct simplex_t {
  int m;
  int n;
  int *var;
  double **a;
  double *b;
  double *x;
  double *c;
  double y;
};


struct node_t{
    int m; 
    int n;
    int k;
    int h;
    double xh;
    double ak;
    double bk;
    double* min;
    double* max;
    int* var;
    double** a; 
    double* b; 
    double* x;
    double* c; 
    double z; 
    struct node_t* next;
};


struct Set {
    struct node_t* first;
};

Set* create_set(){
    Set* s = calloc(1, sizeof(Set));
    s->first = NULL;
    return s;
}

void add_node(Set* s, node_t* p){
    if(s->first == NULL){
        s->first = p;
        s->first->next = NULL;
    }else{
        node_t* n = s->first;
        while(n->next != NULL){
            n = n->next;
        }
        n->next = p;       
    }
}

void remove_node(Set* s, node_t* p){
    if(s->first == NULL){
        return;
    }
    node_t* prev = NULL;
    node_t* n = s->first;
    while(n != p){
        prev = n;
        n = n->next;
        if(n == NULL){
            return;
        }
    }
    prev->next = n->next;
    free(p);

}


node_t* initial_node(int m, int n, double** a, double* b, double* c){
    node_t* p = calloc(1, sizeof(struct node_t)); 
    p->a = calloc(m+1, sizeof(double*)); // matrix
    for(int i = 0; i < m + 1; i++){
        p->a[i] = calloc(n+1, sizeof(double));
    }
    p->b = calloc(m+1, sizeof(double));
    p->c = calloc(n+1, sizeof(double));
    p->x = calloc(n+1, sizeof(double));
    p->min = calloc(n, sizeof(double));
    p->max = calloc(n, sizeof(double));
    p->m = m;
    p->n = n;

    for(int i = 0; i < m; i++){
        p->b[i] = b[i];
        for(int j = 0; j < n; j++){
            p->a[i][j] = a[i][j];
            p->c[j] = c[j];
        }
    }

    for(int i = 0; i < n; i++){
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }
    return p;
}

node_t* extend(node_t* p, int m, int n, double** a, double* b, double* c, int k, double ak, double bk){
    node_t* q = calloc(1, sizeof(struct node_t)); ;
    int i, j;
    q->k = k;
    q->ak = ak;
    q->bk = bk;
    if((ak > 0) && (p->max[k] < INFINITY)){
        q->m = p->m;
    }else if((ak < 0) && (p->min[k] > 0)){
        q->m = p->m;
    }else{
        q->m = p->m +1;
    }
    q->n = p->n;
    q->h = -1;
    q->a = calloc(q->m+1, sizeof(double*)); // matrix
    for(int i = 0; i < q->m + 1; i++){
    
        q->a[i] = calloc(q->n+1, sizeof(double));
    }
    q->b = calloc(q->m+1, sizeof(double));
    q->c = calloc(q->n+1, sizeof(double));
    q->x = calloc(q->n+1, sizeof(double));
    q->min = calloc(n, sizeof(double));
    q->max = calloc(n, sizeof(double));
    for(i = 0; i < n; i++){
        q->min[i] = p->min[i];
        q->max[i] = p->max[i];
    }
    for(i = 0; i < q->m; i ++){
        q->a[i] = a[i];
        q->b[i] = b[i];
    }
    for(i = 0; i < q->n +1; i++){ //q->n +1 ?
        q->c[i] = c[i];
    }
    if(ak > 0){
        if(q->max[k] == INFINITY || bk < q->max[k]){
            q->max[k] = bk;
        } else if(q->min[k] == -INFINITY || -bk > q->min[k]){
            q->min[k] = -bk;
        }
    }
    for(i = m, j = 0; j < n; j++){
        if(q->min[j] > -INFINITY){
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i++;
        }
        if(q->max[j] < INFINITY){
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i++;
        }
    }
    return q;

}

int is_integer(double* xp){
    double x = *xp;
    double r = round(x);
    if(fabs(r-x) < epsilon){
        *xp = r;
        return 1;
    } else {
        return 0;
    }
}

int integer(node_t* p){
    int i;
    for(i = 0; i < p->n; i++){
        if(!is_integer(&p->x[i])){
            return 0;
        }
    }
    return 1;
}

void bound(node_t* p, Set* h, double* zp, double* x){
    if(p->z > *zp){
        *zp = p->z;
        for(int i = 0; i < p->n; i++){
            x[i] = p->x[i];
        }
        node_t* n = h->first;
        while(n!=NULL){
            if(n->z < p->z){
                remove_node(h, n);
            }
            n = n->next;
        }
    }
}
/*
int isfinite(double* x){
    return isfinite(x);
}
*/

int branch(node_t* q, double z){
    double min, max;
    if(q->z < z){
        return 0;
    }
    for(int h = 0; h < q->n; h++){
        if(!is_integer(&q->x[h])){
            if(q->min[h] == -INFINITY){
                min = 0;
            }else{
                min = q->min[h];
            }
            max = q->max[h];
            if(floor(q->x[h])  < min || ceil(q->x[h]) > max){
                continue;
            }
            q->h = h;
            q->xh = q->x[h];
            free(q->a);
            free(q->b);
            free(q->c);
            free(q->x);
            return 1;
        }
    }
    return 0;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y);

void succ(node_t* p, Set* h, int m, int n, double**a, double* b, double*c, int k, double ak, double bk, double* zp, double* x){
    node_t* q = extend(p, m, n, a, b, c, k, ak, bk);
    if(q == NULL){
        return;
    }
    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
    if(isfinite(q->z)){
        if(integer(q)){
            bound(q,h,zp,x);
        }else if(branch(q, *zp)){
            add_node(h, q);
            return;
        }
    }
    free(q);
}

double intopt(double res, int m, int n, double** a, double* b, double* c, double* x){
    node_t* p = initial_node(m,n,a,b,c);
    Set* h = create_set();
    add_node(h, p);
    double z = -INFINITY;
  //  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
    p->z = res;
    if(integer(p) || !isfinite(p->z)){
        z = p->z;
        if(integer(p)){
            for(int i = 0; i < n; i++){
                x[i] = p->x[i];
            }
            
        }
        free(p);
        return z;
    }
    branch(p, p->z);
    node_t* node = h->first;
    while(node!=NULL){
        succ(node, h, m, n, a, b, c, node->h, 1, floor(node->xh), &z, x);
        succ(node, h, m, n, a, b, c, node->h, 1, -ceil(node->xh), &z, x);
        remove_node(h, node);
        node = node->next;
    }
    if(z == - INFINITY){
        return NAN;
    }else{
        return z;
    }
}

//-----------------------------------------------------


int initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var);

int init(simplex_t *s, int m, int n, double **a, double *b, double *c,
         double *x, double y, int *var) {
  int i, k;
  s->m = m;
  s->n = n;
  s->var = var;
  s->a = a;
  s->b = b;
  s->x = x;
  s->c = c;
  s->y = y;

  if (s->var == NULL) {
    s->var = calloc(m + n + 1, sizeof(int));
    for (i = 0; i < m + n; i += 1) {
      s->var[i] = i;
    }
  }
  for (k = 0, i = 1; i < m; i += 1) {
    if (b[i] < b[k]) {
      k = i;
    }
  }
  return k;
}

int select_nonbasic(simplex_t *s) {
  int i;
  for (i = 0; i < s->n; i += 1) {
    if (s->c[i] > epsilon) {
      return i;
    }
  }
  return -1;
}

void pivot(simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;
  int m = s->m;
  int n = s->n;
  int i, j, t;
  t = s->var[col];
  s->var[col] = s->var[n + row];
  s->var[n + row] = t;
  s->y = s->y + c[col] * b[row] / a[row][col];

  for (i = 0; i < n; i += 1) {
    if (i != col) {
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
    }
  }
  c[col] = -c[col] / a[row][col];

  for (i = 0; i < m; i += 1) {
    if (i != row) {
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];
    }
  }
  for (i = 0; i < m; i += 1) {
    if (i != row) {
      for (j = 0; j < n; j += 1) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }
  for (i = 0; i < m; i += 1) {
    if (i != row) {
      a[i][col] = -a[i][col] / a[row][col];
    }
  }
  for (i = 0; i < n; i += 1) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}

void prepare(simplex_t *s, int k) {
  int m = s->m;
  int n = s->n;
  int i;

  for (i = m + n; i > n; i -= 1) {
    s->var[i] = s->var[i - 1];
  }
  s->var[n] = m + n;
  n = n + 1;

  for (i = 0; i < m; i++) {
    s->a[i][n - 1] = -1;
  }
  s->x = calloc(m + n, sizeof(double));
  s->c = calloc(n, sizeof(double));
  s->c[n - 1] = -1;
  s->n = n;
  pivot(s, k, n - 1);
}

double xsimplex(int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h) {
  simplex_t s;
  int i, row, col;
  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    free(s.var);
    return NAN;
  }

  while ((col = select_nonbasic(&s)) >= 0) {
    row = -1;
    for (i = 0; i < m; i += 1) {
      if (a[i][col] > epsilon &&
          (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
        row = i;
      }
    }
    if (row < 0) {
      free(s.var);
      return INFINITY;
    }
    pivot(&s, row, col);
  }

  if (h == 0) {
    for (i = 0; i < n; i += 1) {
      if (s.var[i] < n) {
        x[s.var[i]] = 0;
      }
    }
    for (i = 0; i < m; i += 1) {
      if (s.var[n + 1] < n) {
        x[s.var[n + 1]] = s.b[i];
      }
    }
    free(s.var);
  } else {
    for (i = 0; i < n; i += 1) {
      x[i] = 0;
    }
    for (i = n; i < n + m; i += 1) {
      x[i] = s.b[i - n];
    }
  }
  return s.y;
}

int initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var) {
  int i, j, k;
  double w;
  k = init(s, m, n, a, b, c, x, y, var);
  if (b[k] >= 0) {
    return 1;
  }
  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m + n; i += 1) {
    if (s->var[i] == m + n - 1) {
      if (fabs(s->x[i]) > epsilon) {
        free(s->x);
        free(s->c);
        return 0; // infeasible
      } else {
        break; // This i will be used on the next page
      }
    }
  }

  if (i >= n) {
    // x n+m is basic. find good nonbasic.
    j = 0;
    for (k = 0; k < n; k += 1) {
      if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
        j = k;
      }
    }
    pivot(s, i - n, j);
    i = j;
  }

  if (i < n - 1) {
    // x n+m is nonbasic and not last. swap columns i and n-1
    k = s->var[i];
    s->var[i] = s->var[n - 1];
    s->var[n - 1] = k;
    for (k = 0; k < m; k = k + 1) {
      w = s->a[k][n - 1];
      s->a[k][n - 1] = s->a[k][i];
      s->a[k][i] = w;
    }
  } else {
    // x n+m is nonbasic and last. forget it
  }

  free(s->c);
  s->c = c;
  s->y = y;
  for (k = n - 1; k < n + m - 1; k += 1) {
    s->var[k] = s->var[k + 1];
  }

  n = s->n - 1;
  s->n = s->n - 1;
  double *t;
  t = calloc(n, sizeof(double));

  for (k = 0; k < n; k += 1) {
    for (j = 0; j < n; j += 1) {
      if (k == s->var[j]) {
        t[j] = t[j] + s->c[k];
        goto next_k;
      }
    }
    for (j = 0; j < m; j += 1) {
      if (s->var[n + j] == k) {
        break;
      }
    }
    s->y = s->y + s->c[k] * s->b[j];
    for (i = 0; i < n; i += 1) {
      t[i] = t[i] - s->c[k] * s->a[j][i];
    }
  next_k:;
  }
  for (i = 0; i < n; i += 1) {
    s->c[i] = t[i];
  }
  free(t);
  free(s->x);
  return 1;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

int main(int argc, char** argv)
{
    int m; 
    int n; 
    double** a; 
    double* b;
    double* c; 
    double y = 0;
    double* x; 
    // int* var; 

    // Read rows and cols  
    
    scanf("%d %d8", &m, &n);
    printf("m = %d, n = %d\n", m, n); 

    // Make space for vectors and matrix 
    a = calloc(m+n, sizeof(double*)); // matrix
    for(int i = 0; i < m; i++){
        a[i] = calloc(n+1, sizeof(double));
    }

    b = calloc(m+n, sizeof(double)); //vector 
    c = calloc(n+1, sizeof(double)); // vector 
    //y = calloc(1, sizeof(double));
    x = calloc(n +m+ 1 , sizeof(double)); // might be n 
    // var = calloc(n + m + 1, sizeof(int)); 

    // Scan input 
    for(int i = 0; i < n; i++){
        scanf("%lf", &c[i]);
    } 
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            scanf("%lf", &a[i][j]);
        }
    }
    
    for(int i = 0; i < m; i++){
        scanf("%lf", &b[i]);
    }
    
    double result = simplex(m,n,a,b,c,x,0); 
    double result_intopt = intopt(result, m, n, a, b, c, x);
    printf("Simplex result: %.3lf \n Intopt result: %.3lf", result, result_intopt);


    for(int i = 0; i < m; i++){
        free(a[i]); 
    }
    free(a);  
    free(b);
    free(c);
    free(x); 
    
    return 0; 

}
