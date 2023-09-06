Checking distances within years
================
Amy Bauer,
last edited 2023-09-06

### Step by Step Instructions:

1.  Pick some months that represent the summer (by which I actually mean
    the mosquito period). Say, May-August, but you should choose months
    that make sense given your knowledge of the data and biology.
    Suppose you pick N months to represent the mosquito period. BTW
    these need to be the same months for all locations, so pick a big
    enough window to embrace the activity in all locations.

*-\> see plots at end of document*

``` r
########### Loop through years #########

# define start and end month: period of activity
firstm <- 4
lastm <- 9

# get all unique years
yr.vars <- unique(dat$year)

# data frame to annual values:
all_avg <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(all_avg) <- c("avg_EuDist", "year")
```

2.  For each site and year, you therefore get a length-N vector of
    mosquito abundances for the summer. Normalize each by dividing by
    its Euclidean norm (root mean square of the entries). Let nv(s,y)
    denote this normalized vector for site s and year y.
3.  A measure of intraseasonal synchrony for a year, y, would be to take
    all pairs of *distinct* sites, s1 and s2 and compute the Euclidean
    distance between nv(s1,y) and nv(s2,y). Then average these
    quantities across all pairs of distinct sites.

``` r
for(yr in yr.vars){

  m <- dat %>% 
    filter(year == yr) %>% 
    filter(month >= firstm) %>% 
    filter(month <= lastm) %>% 
    group_by(site) %>% 
    summarize(nv = norm(mpth, type = "2")) %>%  
                   # euclid. norm nv = sqrt(sum(mpth^2))
    select(nv) %>% 
    as.data.frame() %>% 
    distances() %>% # From distances() documentation:
                    #   Euclidean distances
                    #   my_distances1 <- distances(my_data_points)
    distance_columns(1:9) 
  
  # melt matrix into long-format df
  long <- melt(m) %>% rename(V1 = Var1, V2 = Var2) 
  
  # left join to df containing unique site combinations only:
  avg_y <- left_join(s_comb, long) %>% 
    summarize(avg_EuDist = mean(value)) %>% 
    mutate(year = yr)
  
  all_avg <- rbind(all_avg, avg_y)
}

write.csv(all_avg, "output/Cx_tarsalis_annual_avg_eucl_dist.csv", row.names = F)
```

| avg_EuDist | year |
|-----------:|-----:|
|  0.5086066 | 2016 |
|  1.2035905 | 2017 |
|  0.8917539 | 2018 |
|  3.3757705 | 2019 |

*4.To determine if this is significant compared to a null hypothesis of
no phenology whatsoever in mosquito abundances, proceed as follows. For
each s and y, randomize the order of the N entries in the normalized
vector nv(s,y). Then recompute 3 above (whichever of the two metrics in
3 you choose, you have to use the same metric, here, too) to get a
surrogate value. Repeat 10000 times to get a distribution of surrogate
values. The p-value is the fraction of surrogate values smaller (for
metric 1) or larger (form metric 2) than the value computed in 3. Repeat
the whole process for each year.*

------------------------------------------------------------------------

------------------------------------------------------------------------

#### Checking monthly mpth across the 9 sites, 2016-2019

All 9 sites:

![](Cx_tarsalis_yearly_eunorm_eucl_dist_files/figure-gfm/plot%20yearly%20mpth%20by%20site%20over%20month%20I-1.png)<!-- -->

Removing site “OAES” because of the much higher mosquito counts:

![](Cx_tarsalis_yearly_eunorm_eucl_dist_files/figure-gfm/plot%20yearly%20mpth%20by%20site%20over%20month%20II-1.png)<!-- -->
