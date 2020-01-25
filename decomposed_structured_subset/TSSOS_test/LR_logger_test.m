load("model_LR_18.mat")

fname = "LR_18_logger.txt";
model.c = model.C;
run_model_LR(model, "dd", 1, fname)
