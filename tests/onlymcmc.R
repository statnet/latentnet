library(latentnet)

data(sampson)
onlymcmc <- ergmm(samplike ~ latent(d=2, G=3)+rreceiver,tofit="mcmc")

mcmc.diagnostics(onlymcmc)

# Should give an informative error message.
print(try(plot(onlymcmc)))

plot(simulate(onlymcmc))
with(onlymcmc,simulate(model,par=sample[[1]],prior=prior))

heatmap(predict(onlymcmc))