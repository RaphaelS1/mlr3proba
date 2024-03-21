dd = readRDS("~/Downloads/data.rds")
task = dd$task
train_rows = dd$train_rows
test_rows = dd$test_rows

devtools::load_all()

km = lrn("surv.kaplan") #same with cox
km$train(task, row_ids = train_rows)
pred_kaplan = km$predict(task, row_ids = test_rows)

graf_improper = msr("surv.graf", proper = FALSE, id = "graf.improper", p_max = 0.8)
graf_proper   = msr("surv.graf", proper = TRUE,  id = "graf.proper", p_max = 0.8)

km_proper = pred_kaplan$score(graf_proper, task = task, train_set = train_rows)
km_improper = pred_kaplan$score(graf_improper, task = task, train_set = train_rows)

graf_proper$scores
graf_improper$scores
