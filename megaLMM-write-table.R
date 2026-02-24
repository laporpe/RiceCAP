library(rrBLUP)
library(ggplot2)
library(devtools)
library(MegaLMM)
library(dplyr)
library(reshape2) 



# #tutorial Data
# data('yield_data',package='MegaLMM')
# head(yield_data)

#My Data
# rice <- read.csv("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/Phenotype_Data/Filtered_phenotypes_RiceCAPAmp_20231023.csv",header= TRUE,na.strings = c("."))
rice <- read.csv("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/Phenotype_Data/Filtered_phenotypes_RiceCAPAmp_20231023.csv",header= TRUE,na.strings = c("."))
rice[rice=="."] <- NA
head(rice)

rice_filtered <- rice #%>% filter(State == 1 | State ==2)
rice_selected1 <- rice_filtered %>% select(ID,State,REP,TILLNUMMEAN,TOTSEEDWTmean,PANLENGTHPLANTMEAN,MEANHT,SEEDNUMperPANMEAN,
                                           SEEDSAMPLEwtPANmean,TOTBIOMASS,DaysGrainfill,DaysHEAD,DaysMATURITY) %>% filter(DaysHEAD > 0) %>% filter(DaysGrainfill < 150)

rice_selected1$TILLNUMMEAN <- as.numeric(rice_selected1$TILLNUMMEAN)
rice_selected1$TOTSEEDWTmean <- as.numeric(rice_selected1$TOTSEEDWTmean)

#write.csv(rice_selected1, file="/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/Phenotype_Data/BLUPS/rice_selected1.csv")

head(rice_selected1)
str(rice_selected1)

# temperately use genotype mean
rice_geno_mean <- rice_selected1 %>%
  group_by(State, ID) %>%
  summarise_at(c("TILLNUMMEAN","TOTSEEDWTmean","PANLENGTHPLANTMEAN","MEANHT","SEEDNUMperPANMEAN",
                 "SEEDSAMPLEwtPANmean","DaysGrainfill","DaysHEAD","DaysMATURITY","TOTBIOMASS"), mean, na.rm = TRUE) %>% as.data.frame()

head(rice_geno_mean)

# format phenotypic data
library(tibble)
rice_long <- rice_geno_mean %>% 
  melt(id = c("State","ID"))
rice_wide=dcast(rice_long,ID~State+ variable, value.var ="value")
rice_wide=column_to_rownames(.data = rice_wide, var = "ID")
head(rice_long)
head(rice_wide)
dim(rice_wide)
Image(rice_wide)

# My K (is a Kinship matrix... is this correct?)
kinship <- read.table("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/Genotype_Data/Centered_IBS_RiceCAP_Named.txt", header = FALSE, row.names = 1, sep = "\t")
colnames(kinship)=rownames(kinship)
kinship[1:3,1:4]
head(kinship)
dim(kinship)

# prep Y and K
common_names <- intersect(rownames(kinship),rownames(rice_wide))
kinship <- kinship[rownames(kinship) %in% common_names,colnames(kinship) %in% common_names]
rice_wide <- rice_wide %>% filter(rownames(rice_wide) %in% common_names)
rice_wide <- rice_wide[rownames(kinship),]

rice_wide <- rice_wide[,-which(names(rice_wide) %in% c("2_TOTBIOMASS","3_TOTSEEDWTmean","3_TOTBIOMASS","1_TOTBIOMASS"))]

# set up trn/tst - the 3rd column is focal trait

scenario <- c(1:26)

state1 <- c(1:9)
state2 <- c(10:18)
state3 <- c(19:26)

for_state1 <- scenario[!scenario%in%state1]
for_state2 <- scenario[!scenario%in%state2]
for_state3 <- scenario[!scenario%in%state3]

tabler2 <- data.frame(matrix(ncol = 0, nrow = 5))

###############################################################################
###############################################################################
###############################################################################
###############################################################################


traiter <- 1

seed_num <- 1


for(traiter in state3) {
  for(iteration in 1:100) 
  {
    
    trn_dat = tst_dat = rice_wide
    set.seed(seed_num)
    nas = sample(c(rep(TRUE,nrow(rice_wide)/2),rep(FALSE,nrow(rice_wide)/2)))
    
    trait <- colnames(trn_dat)[traiter]
    trn_dat[,traiter][!nas]=NA
    tst_dat[,traiter][nas]=NA
    ###############################################################################
    temp_comp <- c(traiter,for_state3)
    
    
    
    ######trn_dat[which(trait_dat trn_dat[,temp_comp] <- NULL
    trn_dat[,-temp_comp] <- NA
    tst_dat[,-temp_comp] <- NA
    ###############################################################################
    lapply(list(trn_dat, tst_dat, rice_wide), head)
    
    # checking
    sapply(list(rice_wide,kinship), dim)
    sapply(list(rownames(kinship),colnames(kinship),rownames(trn_dat), rownames(tst_dat)), 
           function(x) identical(x,rownames(rice_wide)))
    
    # ------------------------------------------------------------------------------
    # run MegaLMM
    # ------------------------------------------------------------------------------
    
    Y = trn_dat
    K = kinship
    line_data = data.frame(Line = rownames(Y))
    head(line_data)
    
    run_parameters = MegaLMM_control(
      h2_divisions = 20, 
      # Each variance component is allowed to explain between 0% and 100% of the
      # total variation. How many segments should the range [0,100) be divided 
      # into for each random effect?
      burn = 0,  
      # number of burn in samples before saving posterior samples. I set this to 
      # zero and instead run the chain in small chunks, doing the burning manually, a
      # s described below.
      thin = 2,
      # during sampling, we'll save every 2nd sample to the posterior database.
      K = 5 # number of factors. With 19 traits, this is likely way higher than needed.
    )
    
    
    # setAs("matrix", "generalMatrix", function(from) {
    #   Matrix::Matrix(from)
    # })
    
    # Construct the model
    MegaLMM_state = setup_model_MegaLMM(
      Y = Y,  
      # The n x p trait matrix
      formula = ~ 1 + (1|Line),  
      # This is syntax like lme4 for mixed effect models. 
      # We specify a fixed effect of population and a random effect for genotype (Line)
      data = line_data,         
      # the data.frame with information for constructing the model matrices
      relmat = list(Line = K), 
      #######relmat = list(ID = kinship), 
      # A list of covariance matrices to link to the random effects in formula.
      # each grouping variable in formula can be linked to a covariance matrix.
      # If so, every level of the grouping variable must be in the rownames of K.
      # additional rows of K not present in data will still be predicted 
      # (and therefore will use memory and computational time!)
      run_parameters=run_parameters,
      # This list of control parameters created above
      run_ID = "MegaLMM_test"
      # A run identifier. The function will create a folder with this name 
      # and store lots of useful data inside it
    )
    
    # Set priors
    Lambda_prior = list(
      sampler = sample_Lambda_prec_ARD,
      Lambda_df = 3,
      delta_1   = list(shape = 2,  rate = 1),
      delta_2   = list(shape = 3, rate = 1),  # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
      # X = X_Env_Train,
      # X_group = X_Env_groups,
      fit_X = F,  # we start by letting Lambda converge without X, but then turn it on during burnins.
      Lambda_beta_var_shape = 3,
      Lambda_beta_var_rate = 1,
      delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
    )
    
    priors = MegaLMM_priors(
      tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
      tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
      h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
      h2_priors_factors_fun = function(h2s,n) 1, # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
      Lambda_prior = Lambda_prior
    )
    
    # We then assign them to the `MegaLMM_state` object:
    MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
    MegaLMM_state$current_state = NULL
    MegaLMM_state$current_state$delta = NULL
    MegaLMM_state$priors$Lambda_prior$fit_X = FALSE
    MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
    MegaLMM_state$run_parameters$burn = 0
    
    # specify which posterior parameters to save
    MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)
    MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2',
                                                       'tot_Eta_prec','B1','U_F','F',
                                                       'B2_F','Lambda_beta','Lambda_beta_var')
    MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'
    
    MegaLMM_state$Posterior$posteriorFunctions = list(
      Y_Pred  = 'F %*% Lambda + U_R', # total genetic value
      U_Train = '(U_F) %*% Lambda + U_R' # additive genetic value
    )
    # initialize the posterior database
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    
    # Run the MCMC and diagnose convergence
    n_iter = 100
    for(i in 1:10) {
      # if(i > 3) MegaLMM_state$priors$Lambda_prior$fit_X = T
      print(sprintf('Burnin run %d',i))
      # Factor order doesn't "mix" well in the MCMC.
      # We can help it by manually re-ordering from biggest to smallest
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6)
      # clear any previous collected samples because we've re-started the chain 
      MegaLMM_state = clear_Posterior(MegaLMM_state)
      # Draw n_iter new samples, storing the chain
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)
      # # make diagnostic plots
      # traceplot_array(MegaLMM_state$Posterior$Lambda,name = file.path(run_ID,'Lambda.pdf'))
      # # traceplot_array(MegaLMM_state$Posterior$U_Test,name = 'U_Test.pdf',
      # #                 facet_dim = 3)
      # traceplot_array(MegaLMM_state$Posterior$U_Test,name = file.path(run_ID,'U_Test.pdf'))
      # traceplot_array(MegaLMM_state$Posterior$B2_F,facet_dim=3,name = file.path(run_ID,'B2_F.pdf'))
      # traceplot_array(MegaLMM_state$Posterior$Lambda_beta,facet_dim=2,name = file.path(run_ID,'Lambda_beta.pdf'))
      # traceplot_array(MegaLMM_state$Posterior$Lambda_beta_var,facet_dim=2,name = file.path(run_ID,'Lambda_beta_var.pdf'))
      print(sprintf('Completed %d burnin samples', MegaLMM_state$current_state$nrun))
    }
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    
    # Collect posterior samples
    n_iter = 250
    for(i in 1:4) {
      print(sprintf('Sampling run %d',i))
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
      MegaLMM_state = save_posterior_chunk(MegaLMM_state)
      print(MegaLMM_state)
    }
    
    #-------------------------------------------------------------------------------
    ### MegaLMM Predictions
    MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
    #U_Train = get_posterior_mean(MegaLMM_state$Posterior$U_Train)
    #Y_Pred = get_posterior_mean(MegaLMM_state$Posterior$Y_Pred)
    
    U_Train = load_posterior_param(MegaLMM_state,'U_Train')
    Y_Pred = load_posterior_param(MegaLMM_state,'Y_Pred')
    dim(U_Train)
    dim(Y_Pred)
    
    # rrBLUP
    library(rrBLUP)
    gblup=mixed.solve(y=trn_dat[,traiter],K=as.matrix(kinship))
    
    # predictive ability (without estimating gcor)
    temp_Y <- cor(colMeans(Y_Pred[,,traiter]),tst_dat[,traiter],use = "p")
    temp_U <- cor(colMeans(U_Train[,,traiter]),tst_dat[,traiter],use = "p")
    temp_rr <- cor(gblup$u,tst_dat[,traiter],use = "p")
    
    temp_row <- rbind(trait,temp_Y,temp_U,temp_rr,iteration)
    print(temp_row)
    tabler2 <- cbind(tabler2,temp_row)
    
    seed_num = seed_num + 1
    
    # xer <- Y[,traiter][nas]
    # yer = colMeans(U_Train[,,traiter])[nas]
    # 
    # mean_x <- mean(na.omit(xer))
    # mean_y <- mean(na.omit(yer))
    # 
    # # Calculate the intercept for a line with slope = 1 that goes through the mean
    # intercept <- mean_y - mean_x
    # 
    # p <- ggplot(NULL, aes(x=xer, y=yer)) + geom_point() + theme_bw(base_size = 14) + xlab("observed Y") + ylab("expected Y") +
    #   ggtitle(trait)+ annotate("text",  x=Inf, y = Inf, label = paste("  corr: ", round(temp_Y,2)), vjust=2, hjust=1.5) +
    #   geom_abline(slope=1, intercept = intercept, linetype = "dashed", color="Red")
    # 
    # print(p)
    
  }
  
}


write.csv(tabler2, file="/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/Phenotype_Data/BLUPS/MEGALMM_scenario1_100iter_state3.csv")


