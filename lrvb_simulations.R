library(mvtnorm)
library(pracma)
library(ggplot2)

draws_std_normal_1d = function(M=1000){
  xi = rnorm(M, 0, 1) # M-vector
  return(xi)
}

draws_std_normal_2d = function(M=1000){
  xi = rmvnorm(M, mean=rep(0, 2), sigma=diag(2)) # M*2 matrix
  return(xi)
}

get_theta_from_xi_1d = function(xi, mu, sig){
  theta = sapply(xi, function(x){sig*x + mu}) # M-vector
  return(theta)
}

get_theta_from_xi_2d = function(xi, mu, sig){
  theta = t(apply(xi, 1, function(x){sig*x + mu})) # M*K matrix 
  return(theta)
}

KL_mf_1d = function(xi, theta, mu, sig, p, m, Sig, alpha=0){
  Kz = length(m)
  
  get_exact_logpdf_1d = function(x){
    pdf = 1/sqrt(2*pi*Sig)*exp(-1/(2*Sig)*(x-m)^2)
    res = log(sum(p*pdf)) + alpha*x
    # res = log(p[1]*sqrt(1/(2*pi*Sig[1]))*exp(-0.5/Sig[1]*(x-m[1])^2) +
    #   p[2]*sqrt(1/(2*pi*Sig[2]))*exp(-0.5/Sig[2]*(x-m[2])^2)) + alpha*x
    return(res)
  }
  
  kl = - mean(sapply(theta, function(x){get_exact_logpdf_1d(x)})) - log(sig)
  # gradient
  grad_logpdf = sapply(theta, function(x){grad(get_exact_logpdf_1d, x)}) # M-vector
  grad_mu = - mean(grad_logpdf) # K-vector
  grad_sig = - mean(grad_logpdf*xi) - 1/sig
  grad_eta = c(grad_mu, grad_sig)
  # hessian
  hessian_logpdf = sapply(theta, function(x){hessian(get_exact_logpdf_1d, x)}) # M-vector
  grad_mu2 = - mean(hessian_logpdf)
  grad_sig2 = - mean(hessian_logpdf*xi^2) + 1/sig^2
  grad_mu_sig = - mean(hessian_logpdf*xi)
  hessian_eta = matrix(c(grad_mu2, grad_mu_sig, grad_mu_sig, grad_sig2), 2, 2)
  
  KL = list(kl=kl, grad_eta=grad_eta, hessian_eta=hessian_eta)
  return(KL)
}

KL_mf_2d = function(xi, theta, mu, sig, p, m, Sig, alpha=0){
  Kz = length(m) # number of mixtures
  K = length(mu) # dimension
  
  get_exact_logpdf_2d = function(x){
    s1 = Sig[1,1,]
    s2 = Sig[2,2,]
    rho = Sig[1,2,]/sqrt(s1*s2)
    pdf = 1/(2*pi*sqrt(s1*s2*(1-rho^2)))*exp(-0.5/(1-rho^2)*((x[1]-m[,1])^2/s1+
          (x[2]-m[,2])^2/s2-2*rho*(x[1]-m[,1])*(x[2]-m[,2])/sqrt(s1*s2)))
    res = log(sum(p*pdf)) + alpha*x[1]
    return(res)
  }
  
  kl = - mean(apply(theta, 1, function(x){get_exact_logpdf_2d(x)})) - mean(log(sig))
  # gradient
  grad_logpdf = t(apply(theta, 1, function(x){grad(get_exact_logpdf_2d, x)})) # M*2 matrix
  grad_mu = - apply(grad_logpdf, 2, mean) # 2-vector
  grad_sig = - apply(grad_logpdf*xi, 2, mean) - 1/sig # 2-vector
  grad_eta = c(grad_mu, grad_sig) # 4-vector
  # hessian
  hessian_logpdf = t(apply(theta, 1, function(x){hessian(get_exact_logpdf_2d, x)})) # M*4 matrix
  grad_mu2 = - matrix(apply(hessian_logpdf, 2, mean), 2, 2) # 2*2 matirx
  grad_sig2 = - matrix(apply(hessian_logpdf*cbind(xi[,1]^2,xi[,1]*xi[,2],
                      xi[,1]*xi[,2],xi[,2]^2), 2, mean), 2, 2) + diag(1/sig^2)
  grad_mu_sig = - matrix(apply(hessian_logpdf*cbind(xi,xi), 2, mean), 2, 2)
  hessian_eta = matrix(0, 4, 4)
  hessian_eta[1:2,1:2] = grad_mu2
  hessian_eta[3:4,3:4] = grad_sig2
  hessian_eta[1:2,3:4] = hessian_eta[3:4,1:2] = grad_mu_sig
  
  KL = list(kl=kl, grad_eta=grad_eta, hessian_eta=hessian_eta)
  return(KL)
}

KL_opt_1d = function(xi, p, m, Sig, alpha=0, step=1, max_iter=500, tol=1e-6){
  mu = 1
  sig = 2
  eta = c(mu, sig)
  # get initial theta
  theta = get_theta_from_xi_1d(xi, mu, sig)
  # newton update
  kl_value = Inf
  for(j in 1:max_iter){
    kl = KL_mf_1d(xi=xi, theta=theta, mu=mu, sig=sig, p=p, m=m, Sig=Sig, alpha=alpha)
    kl_new = kl$kl
    print(eta)
    print(kl_new)
    J = kl$grad_eta
    H = kl$hessian_eta
    U = solve(H)
    eta_new = eta - step * U %*% J
    if(abs(kl_new - kl_value) < tol){
      break
    }
    eta = eta_new
    mu = eta[1]
    sig = eta[2]
    kl_value = kl_new
    # update theta
    theta = get_theta_from_xi_1d(xi, mu, sig)
  }
  res = list(mean_mf=eta[1], var_mf=eta[2]^2)
  return(res)
}

KL_opt_2d = function(xi, p, m, Sig, alpha=0, step=1, max_iter=500, tol=1e-6){
  mu = c(1, 1)
  sig = c(2, 2)
  eta = c(mu, sig)
  # get initial theta
  theta = get_theta_from_xi_2d(xi, mu, sig)
  # newton update
  kl_value = Inf
  for(j in 1:max_iter){
    kl = KL_mf_2d(xi=xi, theta=theta, mu=mu, sig=sig, p=p, m=m, Sig=Sig, alpha=alpha)
    # print(kl)
    kl_new = kl$kl
    print(eta)
    print(kl_new)
    J = kl$grad_eta
    H = kl$hessian_eta
    U = solve(H)
    eta_new = eta - step * U %*% J
    if(abs(kl_new - kl_value) < tol){
      break
    }
    eta = eta_new
    mu = eta[1:2]
    sig = eta[3:4]
    kl_value = kl_new
    # update theta
    theta = get_theta_from_xi_2d(xi, mu, sig)
  }
  res = list(mean_mf=eta[1:2], var_mf=eta[3:4]^2)
  return(res)
}

LR_var_1d = function(xi, mu, sig, p, m, Sig){
  # mu and sig are estimators of mfvb
  theta = get_theta_from_xi_1d(xi, mu, sig)
  # calculate g_eta
  g_eta = c(1, 0)
  # calculate H_eta2
  kl_at_opt = KL_mf_1d(xi, theta, mu, sig, p, m, Sig)
  H_eta2 = kl_at_opt$hessian_eta
  lr_var = sum(g_eta * solve(H_eta2) %*% g_eta)
  return(lr_var)
}

LR_var_2d = function(xi, mu, sig, p, m, Sig){
  # mu and sig are estimators of mfvb
  theta = get_theta_from_xi_2d(xi, mu, sig)
  # calculate g_eta
  g_eta = c(1, 0, 0, 0)
  # calculate H_eta2
  kl_at_opt_2d = KL_mf_2d(xi, theta, mu, sig, p, m, Sig)
  H_eta = kl_at_opt_2d$hessian_eta
  lr_var = sum(g_eta * solve(H_eta) %*% g_eta)
  return(lr_var)
}

lap_mean_var_1d = function(p, m, Sig, alpha=0){
  fun1 = function(x){
    pdf = 1/sqrt(2*pi*Sig)*exp(-1/(2*Sig)*(x-m)^2)
    res = log(sum(p*pdf)) + alpha*x
    # res = log(p[1]*sqrt(1/(2*pi*Sig[1]))*exp(-0.5/Sig[1]*(x-m[1])^2) +
    #   p[2]*sqrt(1/(2*pi*Sig[2]))*exp(-0.5/Sig[2]*(x-m[2])^2)) + alpha*x
    return(res)
  }
  lap = optimize(fun1, c(-2, 2), maximum=T)
  mean_lap = lap$maximum
  # calculate variance
  hess_at_max = hessian(fun1, mean_lap)
  var_lap = - 1 / hess_at_max
  res = list(mean_lap=mean_lap, var_lap=var_lap)
  return(res)
}

lap_mean_var_2d = function(p, m, Sig, alpha=0){
  fun2 = function(x){
    s1 = Sig[1,1,]
    s2 = Sig[2,2,]
    rho = Sig[1,2,]/sqrt(s1*s2)
    pdf = 1/(2*pi*sqrt(s1*s2*(1-rho^2)))*exp(-0.5/(1-rho^2)*((x[1]-m[,1])^2/s1+
        (x[2]-m[,2])^2/s2-2*rho*(x[1]-m[,1])*(x[2]-m[,2])/sqrt(s1*s2)))
    res = log(sum(p*pdf)) + alpha*x[1]
    return(-res) # make is a minimization problem
  }
  lap = optim(c(0, 0), fun2, hessian=T)
  mean_lap = lap$par
  # calculate variance
  hess_at_max = lap$hessian
  var_lap = solve(hess_at_max)
  res = list(mean_lap=mean_lap, var_lap=var_lap)
  return(res)
}

draws_mixture_1d = function(p, m, Sig, M=1000){
  # Kz-component mixture model
  I = rmultinom(M, size=1, prob=p)
  d = apply(I, 2, function(x){
    k = which(x==1)
    rnorm(1, mean=m[k], sd=sqrt(Sig[k]))
  })
  return(d)
}

draws_mixture_2d = function(p, m, Sig, M=1000){
  I = rmultinom(M, size=1, prob=p)
  d = apply(I, 2, function(x){
    k = which(x==1)
    rmvnorm(1, mean=m[k,], sigma=Sig[,,k])
  }) # 2*M matrix
  return(t(d)) # M*2 matrix
}

mc_mean_var_1d = function(p, m, Sig, M=1000){
  draws = draws_mixture_1d(p, m, Sig, M)
  mean_mc = mean(draws)
  var_mc = mean(draws^2) - mean_mc^2
  res = list(mean_mc=mean_mc, var_mc=var_mc)
  return(res)
}

mc_mean_var_2d = function(p, m, Sig, M=1000){
  draws = draws_mixture_2d(p, m, Sig, M) # M*2 matrix
  mean_mc = apply(draws, 2, mean) # 2-vector
  var_mc = t(draws)%*%draws/M - mean_mc%*%t(mean_mc) # 2*2 matrix
  maginal_var_mc = c(var_mc[1,1], var_mc[2,2]) # 2-vector
  res = list(mean_mc=mean_mc, var_mc=var_mc)
  return(res)
}

##### univariate skewed distribution #####
# 2-component
set.seed(1)
p = c(0.5, 0.5)
m = c(0, 2.5) # true mean
Sig = c(2/3, 2) # true variance
alpha = 0
M = 1000

xi = draws_std_normal_1d(M)

# compute mean-field mean and variance
res_mf = KL_opt_1d(xi, p, m, Sig, alpha=0, step=0.1)
mean_mf = res_mf$mean_mf
var_mf = res_mf$var_mf
sig_mf = sqrt(var_mf)

# compute exact mean and variance using MC
res_mc = mc_mean_var_1d(p, m, Sig, M)
mean_mc = res_mc$mean_mc
var_mc = res_mc$var_mc

# compute laplace mean and variance
res_lap = lap_mean_var_1d(p, m, Sig, alpha)
mean_lap = res_lap$mean_lap
var_lap = res_lap$var_lap

# compute lrvb variance
var_lr = LR_var_1d(xi, mean_mf, sig_mf, p, m, Sig)

# result
res_uni_skew = list(Metric=c('mean','var'),Exact=c(mean_mc,var_mc),
                    LRVB=c(mean_mf,var_lr),MFVB=c(mean_mf,var_mf),
                    Laplace=c(mean_lap,var_lap))
res_uni_skew

##### plot #####
# exact posterior and approximations
x1 = seq(-3, 6, length.out=1000)
exact_pdf_uni = function(x, alpha){
  pdf = 1/sqrt(2*pi*Sig)*exp(-1/(2*Sig)*(x-m)^2)
  res = sum(p*pdf) * exp(alpha*x)
  return(res)
}
alpha = 0
p_exact = sapply(x1, exact_pdf_uni, alpha)
p_exact = p_exact/(sum(p_exact)*9/1000) # normalize
p_mfvb = dnorm(x1, mean=mean_mf, sd=sqrt(var_mf))
p_mfvb = p_mfvb/(sum(p_mfvb)*9/1000)
p_lap = dnorm(x1, mean=mean_lap, sd=sqrt(var_lap))
p_lap = p_lap/(sum(p_lap)*9/1000)
plot_uni_skew = data.frame(x=x1, Exact=p_exact, MFVB=p_mfvb, Laplace=p_lap)
png('uni_skew_pdf.png')
ggplot(data=plot_uni_skew) +
  geom_line(aes(x=x1, y=Exact, color='Exact'), size=1.5) +
  geom_line(aes(x=x1, y=MFVB, color='MFVB'), size=1.5) +
  geom_line(aes(x=x1, y=Laplace, color='Laplace'), size=1.5) +
  geom_vline(aes(xintercept=mean_mc, color='Exact'), size=1.5) +
  geom_vline(aes(xintercept=mean_mf, color='MFVB'), size=1.5) +
  geom_vline(aes(xintercept=mean_lap, color='Laplace'), size=1.5) +
  ggtitle('Exact posterior and approximations') +
  xlab(expression(theta[1])) +
  ylab('Density') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# effect of tilting on exact posterior
alpha = 0.1
p_exact_01 = sapply(x1, exact_pdf_uni, alpha)
p_exact_01 = p_exact_01/(sum(p_exact_01)*9/1000)
tilt_pdf = data.frame(x=x1, alpha0=p_exact, alpha1=p_exact_01)
mean_exact = sum(x1*p_exact)*9/1000
mean_exact_01 = sum(x1*p_exact_01)*9/1000
png('uni_skew_tilt_pdf.png')
ggplot(data=tilt_pdf) +
  geom_line(aes(x=x1, y=alpha0, color='alpha0', lty='alpha0'), size=1.5) +
  geom_line(aes(x=x1, y=alpha1, color='alpha1', lty='alpha1'), size=1.5) +
  geom_vline(aes(xintercept=mean_exact, color='alpha0', lty='alpha0'), size=1.5) +
  geom_vline(aes(xintercept=mean_exact_01, color='alpha1', lty='alpha1'), size=1.5) +
  ggtitle('Effect of tilting on exact posterior') +
  xlab(expression(theta[1])) +
  ylab('Density') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# posterior means vs amount of tilting
res_mf_01 = KL_opt_1d(xi, p, m, Sig, alpha=0.1, step=0.1)
mean_mf_01 = res_mf_01$mean_mf
res_lap_01 = lap_mean_var_1d(p, m, Sig, alpha=0.1)
mean_lap_01 = res_lap_01$mean_lap
slope_exact = (mean_exact_01 - mean_exact)/0.1
slope_mf = (mean_mf_01 - mean_mf)/0.1
slope_lap = (mean_lap_01 - mean_lap)/0.1
plot_mean = data.frame(method=c('Exact', 'MFVB', 'Laplace'), 
                       slope=c(slope_exact, slope_mf, slope_lap))
png('uni_skew_mean_tilt.png')
ggplot(data=plot_mean) +
  geom_abline(aes(intercept=0, slope=slope, color=method, lty=method), size=1.5) +
  xlim(0, 0.1) +
  ylim(0, 0.35) +
  ggtitle('Posterior means vs amount of tilting') +
  xlab(expression(alpha)) +
  ylab('Posterior mean') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()


##### univariate over-dispersed distribution #####
set.seed(2)
p = c(1/3, 1/3, 1/3)
m = c(-2, 0, 2)
Sig = c(2, 2/3, 2)
alpha = 0
M = 1000

xi = draws_std_normal_1d(M)

# compute mean-field mean and variance
res_mf = KL_opt_1d(xi, p, m, Sig, alpha, step=0.1)
mean_mf = res_mf$mean_mf
var_mf = res_mf$var_mf
sig_mf = sqrt(var_mf)

# compute exact mean and variance using MC
res_mc = mc_mean_var_1d(p, m, Sig, M)
mean_mc = res_mc$mean_mc
var_mc = res_mc$var_mc

# compute laplace mean and variance
res_lap = lap_mean_var_1d(p, m, Sig, alpha)
mean_lap = res_lap$mean_lap
var_lap = res_lap$var_lap

# compute lrvb variance
var_lr = LR_var_1d(xi, mean_mf, sig_mf, p, m, Sig)

# result
res_uni_over = list(Metric=c('mean','var'),Exact=c(mean_mc,var_mc),
                    LRVB=c(mean_mf,var_lr),MFVB=c(mean_mf,var_mf),
                    Laplace=c(mean_lap,var_lap))
res_uni_over

##### plot #####
# exact posterior and approximations
x1 = seq(-6, 6, length.out=1000)
alpha = 0
p_exact = sapply(x1, exact_pdf_uni, alpha)
p_exact = p_exact/(sum(p_exact)*12/1000) # normalize
p_mfvb = dnorm(x1, mean=mean_mf, sd=sqrt(var_mf))
p_mfvb = p_mfvb/(sum(p_mfvb)*12/1000)
p_lap = dnorm(x1, mean=mean_lap, sd=sqrt(var_lap))
p_lap = p_lap/(sum(p_lap)*12/1000)
plot_uni_over = data.frame(x=x1, Exact=p_exact, MFVB=p_mfvb, Laplace=p_lap)
png('uni_over_pdf.png')
ggplot(data=plot_uni_over) +
  geom_line(aes(x=x1, y=Exact, color='Exact'), size=1.5) +
  geom_line(aes(x=x1, y=MFVB, color='MFVB'), size=1.5) +
  geom_line(aes(x=x1, y=Laplace, color='Laplace'), size=1.5) +
  geom_vline(aes(xintercept=mean_mc, color='Exact'), size=1.5) +
  geom_vline(aes(xintercept=mean_mf, color='MFVB'), size=1.5) +
  geom_vline(aes(xintercept=mean_lap, color='Laplace'), size=1.5) +
  ggtitle('Exact posterior and approximations') +
  xlab(expression(theta[1])) +
  ylab('Density') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# effect of tilting on exact posterior
alpha = 0.1
p_exact_01 = sapply(x1, exact_pdf_uni, alpha)
p_exact_01 = p_exact_01/(sum(p_exact_01)*12/1000)
tilt_pdf = data.frame(x=x1, alpha0=p_exact, alpha1=p_exact_01)
mean_exact = sum(x1*p_exact)*12/1000
mean_exact_01 = sum(x1*p_exact_01)*12/1000
png('uni_over_tilt_pdf.png')
ggplot(data=tilt_pdf) +
  geom_line(aes(x=x1, y=alpha0, color='alpha0', lty='alpha0'), size=1.5) +
  geom_line(aes(x=x1, y=alpha1, color='alpha1', lty='alpha1'), size=1.5) +
  geom_vline(aes(xintercept=mean_exact, color='alpha0', lty='alpha0'), size=1.5) +
  geom_vline(aes(xintercept=mean_exact_01, color='alpha1', lty='alpha1'), size=1.5) +
  ggtitle('Effect of tilting on exact posterior') +
  xlab(expression(theta[1])) +
  ylab('Density') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# posterior means vs amount of tilting
res_mf_01 = KL_opt_1d(xi, p, m, Sig, alpha=0.1, step=0.1)
mean_mf_01 = res_mf_01$mean_mf
res_lap_01 = lap_mean_var_1d(p, m, Sig, alpha=0.1)
mean_lap_01 = res_lap_01$mean_lap
slope_exact = (mean_exact_01 - mean_exact)/0.1
slope_mf = (mean_mf_01 - mean_mf)/0.1
slope_lap = (mean_lap_01 - mean_lap)/0.1
plot_mean = data.frame(method=c('Exact', 'MFVB', 'Laplace'), 
                       slope=c(slope_exact, slope_mf, slope_lap))
png('uni_over_mean_tilt.png')
ggplot(data=plot_mean) +
  geom_abline(aes(intercept=0, slope=slope, color=method, lty=method), size=1.5) +
  xlim(0, 0.1) +
  ylim(0, 0.45) +
  ggtitle('Posterior means vs amount of tilting') +
  xlab(expression(alpha)) +
  ylab('Posterior mean') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()


##### bivariate over-dispersed distribution #####
set.seed(3)
M = 1000
w_diff = 0.15
p = c(1/3-0.5*w_diff, 1/3+w_diff, 1/3-0.5*w_diff)
delta = 1.5
m = matrix(1,3,2)
m[1,] = -delta*m[1,]
m[2,] = 0
m[3,] = delta*m[3,]
sd_narrowness = 0.4
sd_scale = 0.6
get_rotation_matrix = function(x){
  return(matrix(c(cos(x),sin(x),-sin(x),cos(x)), 2, 2))
}
rotate_mat = function(mat, x){
  rot_mat = get_rotation_matrix(x)
  return(rot_mat%*%mat%*%solve(rot_mat))
}
Sig0 = diag(c(1,sd_narrowness^2))*(sd_scale*delta)^2
Sig = rotate_mat(Sig0, pi/4)
Sig = replicate(3, Sig, simplify='array')
alpha = 0

xi = draws_std_normal_2d(M)

# compute mean-field mean and variance
res_mf = KL_opt_2d(xi, p, m, Sig, alpha, step=0.1)
mean_mf = res_mf$mean_mf
var_mf = res_mf$var_mf
sig_mf = sqrt(var_mf)

# compute exact mean and variance using MC
res_mc = mc_mean_var_2d(p, m, Sig, M)
mean_mc = res_mc$mean_mc
var_mc = res_mc$var_mc

# compute laplace mean and variance
res_lap = lap_mean_var_2d(p, m, Sig, alpha)
mean_lap = res_lap$mean_lap
var_lap = res_lap$var_lap

# compute lrvb variance
var_lr = LR_var_2d(xi, mean_mf, sig_mf, p, m, Sig)

# result
res_bi_over = list(Metric=c('mean','var'),Exact=c(mean_mc[1],var_mc[1,1]),
                    LRVB=c(mean_mf[1],var_lr),MFVB=c(mean_mf[1],var_mf[1]),
                    Laplace=c(mean_lap[1],var_lap[1,1]))
res_bi_over

##### plot #####
# exact posterior and approximations
x1 = seq(-4, 4, length.out=100)
x1 = as.matrix(expand.grid(x=x1, y=x1))
exact_pdf_bi_over = function(x, alpha){
  s1 = Sig[1,1,]
  s2 = Sig[2,2,]
  rho = Sig[1,2,]/sqrt(s1*s2)
  pdf = 1/(2*pi*sqrt(s1*s2*(1-rho^2)))*exp(-0.5/(1-rho^2)*((x[1]-m[,1])^2/s1+
                                                             (x[2]-m[,2])^2/s2-2*rho*(x[1]-m[,1])*(x[2]-m[,2])/sqrt(s1*s2)))
  res = sum(p*pdf)*exp(alpha*x[1])
  return(res)
}
alpha = 0
p_exact = apply(x1, 1, exact_pdf_bi_over, alpha)
p_exact = p_exact/(sum(p_exact)*16/10000) # normalize
p_mfvb = dmvnorm(x1, mean=mean_mf, sigma=diag(var_mf))
p_mfvb = p_mfvb/(sum(p_mfvb)*16/10000)
p_lap = dmvnorm(x1, mean=mean_lap, sigma=var_lap)
p_lap = p_lap/(sum(p_lap)*16/10000)
plot_bi_over = data.frame(x=x1, Exact=p_exact, MFVB=p_mfvb, Laplace=p_lap)
png('bi_over_contour.png')
ggplot(data=plot_bi_over) +
  geom_contour(aes(x=x1[,1], y=x1[,2], z=Exact, color='Exact'), size=1.5) +
  geom_contour(aes(x=x1[,1], y=x1[,2], z=MFVB, color='MFVB'), size=1.5) +
  geom_contour(aes(x=x1[,1], y=x1[,2], z=Laplace, color='Laplace'), size=1.5) +
  ggtitle('Exact posterior and approximations') +
  xlab(expression(theta[1])) +
  ylab(expression(theta[2])) +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# marginal distribution and approximations
marginal_bi_over = function(x, alpha){
  pdf = 1/sqrt(2*pi*Sig[1,1,])*exp(-1/(2*Sig[1,1,])*(x-m[,1])^2)
  res = sum(p*pdf) * exp(alpha*x)
  return(res)
}
x1 = seq(-4, 4, length.out=1000)
p_exact = sapply(x1, marginal_bi_over, alpha)
p_exact = p_exact/(sum(p_exact)*8/1000) # normalize
p_mfvb = dnorm(x1, mean=mean_mf[1], sd=sqrt(var_mf[1]))
p_mfvb = p_mfvb/(sum(p_mfvb)*8/1000)
p_lap = dnorm(x1, mean=mean_lap[1], sd=sqrt(var_lap[1,1]))
p_lap = p_lap/(sum(p_lap)*8/1000)
plot_bi_over = data.frame(x=x1, Exact=p_exact, MFVB=p_mfvb, Laplace=p_lap)
png('bi_over_marginal.png')
ggplot(data=plot_bi_over) +
  geom_line(aes(x=x1, y=Exact, color='Exact'), size=1.5) +
  geom_line(aes(x=x1, y=MFVB, color='MFVB'), size=1.5) +
  geom_line(aes(x=x1, y=Laplace, color='Laplace'), size=1.5) +
  geom_vline(aes(xintercept=mean_mc[1], color='Exact'), size=1.5) +
  geom_vline(aes(xintercept=mean_mf[1], color='MFVB'), size=1.5) +
  geom_vline(aes(xintercept=mean_lap[1], color='Laplace'), size=1.5) +
  ggtitle('Marginal distribution and approximations') +
  xlab(expression(theta[1])) +
  ylab('Density') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()

# posterior means vs amount of tilting
p_exact_01 = sapply(x1, marginal_bi_over, alpha=0.1)
p_exact_01 = p_exact_01/(sum(p_exact_01)*8/1000)
mean_exact = sum(x1*p_exact)*8/1000
mean_exact_01 = sum(x1*p_exact_01)*8/1000
res_mf_01 = KL_opt_2d(xi, p, m, Sig, alpha=0.1, step=0.1)
mean_mf_01 = res_mf_01$mean_mf
res_lap_01 = lap_mean_var_2d(p, m, Sig, alpha=0.1)
mean_lap_01 = res_lap_01$mean_lap
slope_exact = (mean_exact_01[1] - mean_exact[1])/0.1
slope_mf = (mean_mf_01[1] - mean_mf[1])/0.1
slope_lap = (mean_lap_01[1] - mean_lap[1])/0.1
plot_mean = data.frame(method=c('Exact', 'MFVB', 'Laplace'), 
                       slope=c(slope_exact, slope_mf, slope_lap))
png('bi_over_mean_tilt.png')
ggplot(data=plot_mean) +
  geom_abline(aes(intercept=0, slope=slope, color=method, lty=method), size=1.5) +
  xlim(0, 0.01) +
  ylim(0, 0.018) +
  ggtitle('Posterior means vs amount of tilting') +
  xlab(expression(alpha)) +
  ylab('Posterior mean') +
  theme(plot.title=element_text(hjust=0.5))
dev.off()


