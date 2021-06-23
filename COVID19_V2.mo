model COVID19_V2 "SEIR Compartmental Model"
  parameter Real u=0.2 "social distancing term (u=0 (none), u=1 (total isolation))";
  parameter Real t_incubation = 5.1*24*3600 "Incubation period";
  parameter Real t_infective = 3.3*24*3600 "Infective period";
  parameter Real R0 = 2.4 "Basic reproduction ratio";
  parameter Real N = 100000 "Total population";
  // initial number of infected and recovered individuals
  parameter Real e_initial = 1/N;
  parameter Real i_initial = 0.00;
  parameter Real r_initial = 0.00;
  parameter Real s_initial = 1 - e_initial - i_initial - r_initial;
  parameter Real alpha = 1/t_incubation;
  parameter Real gamma = 1/t_infective;
  parameter Real beta = R0*gamma;
  Real s(start=s_initial)"Number susceptible";
  Real e(start=e_initial)"Number exposed";
  Real i(start=i_initial)"Number infectious";
  Real r(start=r_initial)"Number recovered";
equation
  der(s)=-(1-u)*beta*s*i;
  der(e)=(1-u)*beta*s*i-alpha*e;
  der(i)=alpha*e-gamma*i;
  der(r)=gamma*i;
end COVID19_V2;
