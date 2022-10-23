packages <- c("generics", "demography", "forecast")
## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#### load one-way median polish function
source("one_way_median_polish.R")
###read data
Australia_data =  demography::extract.ages(read.demogdata("Australia_mort.txt", "Australia_expo.txt", type = "mortality", label = "AUS"),
                               ages = 0:100)
Australia_female = t(log(Australia_data$rate$female, base = 10))
year = Australia_data$year
n_year = length(year)
age = Australia_data$age
n_age = length(age)

# create a `not-so-useful` list, unsure if it is correct

part_list = list()
for(ik in 1:n_year) {
  part_list[[ik]] = ik
}

# one-way functional data (matrix)
Australia_female_median_polish = one_way_median_polish(Y = Australia_female, partition_index = part_list)
grand_effect <- Australia_female_median_polish$grand_effect
row_effect <- Australia_female_median_polish$row_effect
# residuals

Australia_female_resi <- sapply(1:n_year, function(k) {
  value <- Australia_female[k, ] - grand_effect-row_effect[k, ]
  return (value)
})


# two populations (female and male)
source("one_way_median_polish.R")
Australia_both = array(NA, dim = c(n_year, n_age, 2))
Australia_both[, , 1] = t(log(Australia_data$rate$female, base = 10))
Australia_both[, , 2] = t(log(Australia_data$rate$male, base = 10))


Australia_both_median_polish <- one_way_median_polish(Y = Australia_both, 
                                                      partition_index = part_list)
grand_effect <- Australia_both_median_polish$grand_effect
row_effect <- Australia_both_median_polish$row_effect

# below produces an error 
Australia_both_resi <- Australia_both
for (k in 1:n_year) {
  Australia_both_resi[k, , ] <- Australia_both[k, , ] - grand_effect - row_effect[[k]]
}

