Simulate UMI counts
===================

PROSSTT simulates UMI counts using a negative binomial distribution where the variance of the expression of each gene :math:`g` depends on its average expression :math:`\mu_g`:

:math:`\sigma_g^2 = \alpha_g \mu^2 + \beta_g \mu_g`

This relationship is preserved through pseudotime; as average expression changes with time, so does the variance, always obeying the same relationship.

:math:`\sigma_g^2(t) = \alpha_g \mu^2(t) + \beta_g \mu_g(t)`

A negative binomial is the distribution of the number of successes in a sequence of i.i.d. Bernoulli trials before a specified number of failures occurs. The negative binomial can be parametrized by its mean and variance or by a pair :math:`p \in (0, 1), r > 0`, where :math:`p` is the success probability in each Bernoulli trial and :math:`r` the number of failures. While the negative binomial is originally a discrete probability distribution, it can easily be extended into a continuous one, preserving most of its attributes.

Here we use the implementation of the negative binomial distribution by the `scipy` package, after translating the mean and variance of each distribution to the :math:`p, r` equivalents.

The gene-specific parameters :math:`\alpha_g, \beta_g` are sampled from ranges found in real data. Users can set :math:`\alpha_g` to 0 and :math:`\beta_g` to 1 to have genes with Poisson distributions, or only set :math:`\alpha_g` to 0 to have genes with scaled Poissonian noise.