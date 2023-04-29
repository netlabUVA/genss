gen_rinit_lin <- function() {
  return(pomp::Csnippet("
          state1 = rnorm(0,1);
          state2 = rnorm(0,1);
        "))
}

gen_rprocess_lin <- function() {
  return(pomp::Csnippet("
          double d = (a1_1*a2_1+a1_2*a2_2)/(1-a1_1*a2_2 -a1_2*a2_1);
          double state1_var = 1 -(a1_1*a1_1 +a1_1*a1_2*d*2 + a1_2*a1_2);
          double state2_var = 1 -(a2_1*a2_1 +a2_1*a2_2*d*2 + a2_2*a2_2);
          
          if(state1_var < 0){
          state1_var = .5;
          }
          if(state2_var < 0){
          state2_var = .5;
          }
          double tstate1 = rnorm(a1_1*state1 + a1_2*state2, sqrt(state1_var));
          double tstate2 = rnorm(a2_1*state1 + a2_2*state2, sqrt(state2_var));
          state1 = tstate1;
          state2 = tstate2;
        "))
}
