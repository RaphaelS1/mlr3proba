t = tsk("rats")
s = partition(t)
p = lrn("surv.coxph")$train(t, s$train)$predict(t, s$test)
m_proper = msr("surv.graf", proper = TRUE)
m_improper = msr("surv.graf", proper = FALSE)
p$score(m_proper)
p$score(m_improper)

mean(abs(m_proper$scores - m_improper$scores))
mean((m_proper$scores - m_improper$scores)^2)

plot(m_proper$scores, m_improper$scores)
abline(a = 0, b = 1)
