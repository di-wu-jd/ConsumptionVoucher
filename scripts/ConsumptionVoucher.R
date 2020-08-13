# Title     : Scripts to process experiment data on Consumer Voucher
# Updated on: 8/8/20

library(data.table)
library(WRS)

# Loading data
data <- fread('data/consumption_vouchers_jd.csv')

# data inspection
# total number of rows
data[, .N]

# data sample
head(data)

# Basic statistics
data <- data[, converted:=0]
data <- data[gmv>0, converted:=1]
data <- data[, converted_total:=0]
data <- data[gmv_total>0, converted_total:=1]
# user counts, means, trimmed means, conversion rates and take up rates
data[, list(user_count=.N,
            mean_gmv=mean(gmv),
            trimmed_mean_gmv=mean(gmv, trim=0.001),
            conversion_rate=sum(converted)/.N,
            take_up_rate=sum(used_coupon)/.N,
            long_term_mean_gmv=mean(gmv_total),
            long_term_trimmed_mean_gmv=mean(gmv_total, trim=0.001),
            long_term_conversion_rate=sum(converted_total)/.N),
       by=group][order(group)]

# GMV percentiles
data[, list(user_count=.N, mean_gmv=mean(gmv), gmv_se=sd(gmv),
            percentile_99=quantile(gmv, c(.99)),
            percentile_991=quantile(gmv, c(.991)),
            percentile_992=quantile(gmv, c(.992)),
            percentile_993=quantile(gmv, c(.993)),
            percentile_994=quantile(gmv, c(.994)),
            percentile_995=quantile(gmv, c(.995)),
            percentile_996=quantile(gmv, c(.996)),
            percentile_997=quantile(gmv, c(.997)),
            percentile_998=quantile(gmv, c(.998)),
            percentile_999=quantile(gmv, c(.999)),
            percentile_9995=quantile(gmv, c(.9995))), by=group][order(group)]

# Bootstrapping

# preprocessing data
x <- fac2list(data$gmv, data$group)
# set up comparison
con <- matrix(0, nrow=16, ncol=15)
con[1,] <- -1
for (i in seq(1, 15)) {
  con[i+1, i] <- 1
}
# run bootstrapping with trimmed means (trim=0.001)
linconpba(x=x, con=con, nboot=2000, ni=400, alpha=0.05, est=mean, bhop=FALSE, SEED=TRUE, trim=0.001)

# run bootstrap to calculate confidence interval for returns on trimmed data.
# Example: Compare group 5 v.s. control
boot_return(c=data[group == 0, .(outcome=gmv, cost=coupon_cost*used_coupon)],
            t=data[group == 5, .(outcome=gmv, cost=coupon_cost*used_coupon)],
            tr=0.001, nboot=2000, ni=400, SEED=TRUE, alpha=0.05)

