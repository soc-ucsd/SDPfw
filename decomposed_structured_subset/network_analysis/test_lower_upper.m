%load('sea_star_Hinf0_small.mat', 'model_lower', 'model_upper')
load('sea_star_Hinf0_verylarge.mat', 'model_lower', 'model_upper')

prob_lower = convert_sedumi2mosek(model_lower.A, model_lower.b, model_lower.c, model_lower.K);
prob_upper = convert_sedumi2mosek(model_upper.A, model_upper.b, model_upper.c, model_upper.K);

[r_lower, res_lower] = mosekopt('minimize', prob_lower);

cost_lower = res_lower.sol.itr.pobjval;

[r_upper, res_upper] = mosekopt('minimize', prob_upper);

cost_upper = res_upper.sol.itr.pobjval;