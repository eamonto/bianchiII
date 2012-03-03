///// STORE THE CONDITIONS BEFORE THE INTEGRATION
int store_levels_rk45(void);

///// FINAL VALUES AFTER INTEGRATION
int step_rk45(int i);

///// EVOLUTION WITH THE ADAPTIVE METHOD
int evolution_rk45(int i);

///// RELATIVE ERROR BETWEEN THE RK4 AND RK5 METHODS
int error_rk45(void);
