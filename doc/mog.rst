.. _mog:

Mixture-of-Gaussians EM Algorithm
==================================

This software uses KlustaKwik_ to do the clustering on subsets of channels. I originally used a different approach to clustering, where all channels were clustered at once using a modified EM algorithm. When I switched back to using the standard mixture-of-Gaussians EM algorithm on subsets of channels, I switched to KlustaKwik, since it is faster and better tested than my clustering program (which can be found in ``core/CEM_adjacency.py``.)

I wrote up these notes as part of the explanation of that clustering method.  Now this section is possibly unnecessary, but I include it just in case someone is interested.

EM Algorithm
--------------------

Clustering can be cast as a statistical inference problem like linear regression, in which we find the best model to fit our data. Here, the model is a *mixture model*, where there are several different mixture components with different parameters producing data points. The data-generating process in a mixture-of-Gaussians model is as follows: 

   Repeat this N times independently, where N is the number of data points: 

   *  Randomly select a mixture component m with probability  :math:`\pi(m)`
   *  Select a point from the probability distribution of this component, using probability distribution :math:`P(X|m)=\mathcal{N}(X|\mu_{m},\Lambda_{m})`, where :math:`\mu` and :math:`\Lambda` are means and covariances.

We assume that our dataset, i.e., our collection of feature vectors, was generated using this process. Our goal is to fit the parameters :math:`\mu,\Lambda,\pi` to the data. In particular, we look for the *maximum likelihood* values of the parameters, which are the parameter values that maximize P(Data | Parameters). From now on, we will use X as shorthand for "data", and :math:`\Theta` as shorthand for "parameters", i.e., :math:`\{\mu,\Lambda,\pi\}`. Note that :math:`P(X | \Theta)` is a probability *density*, not a finite probability. 

Suppose Z is the latent variable that says which observed vector comes from which cluster, called the *responsibilities*. Then the *likelihood* L is as follows:

.. math:: L = P(X|\mu,\Lambda,\pi)=\sum_{Z} P(X|Z,\mu,\Lambda)P(Z|\pi)

We iterate the following two steps to maximize L over the parameters: 

   *  E Step: update Z, i.e., or more precisely, :math:`P(Z | X,\mu,\Lambda,\pi)`.
   *  M Step: update the maximum likelihood values :math:`\pi, \mu, \Lambda` for each cluster.
   
We compute likelihood score :math:`L = p(X | \Theta)` after each step. L is guaranteed to increase after each EM step pair. This is because the likelihood can be written as

.. math::   \log P(X|\Theta) = E_{Q(Z)} P(X|\Theta,Z) + KL(Q(Z) \| P(Z|X,\Theta))
   
KL is the Kullback-Liebler divergence, a nonnegative distance between two probability distributions. In words, the above equation reads

.. math:: 
   :nowrap:
   
   \begin{align*}
      \log P(X|\Theta) = (&\text{log of the expected value of P(X) given some wrong distribution Q(Z) over the latent variable Z}) \\
      + (&\text{the mismatch between Q(Z) and the correct probability distribution } P(Z|X,\Theta)) 
   \end{align*}

The likelihood score is the first term on the right-hand side, :math:`E_{Q(Z)} P(X|\Theta,Z)`. In the M step, we strictly increase the first term :math:`E_{Q(Z)} P(X|\Theta,Z)` while keeping the second term :math:`KL(Q(Z) \| P(Z|X,\Theta))` constant. On the E step, we drive the second term to zero, but the sum :math:`\log P(X|\Theta)` remains constant, so the first term, the likelihood score, must decrease.

This is based on chapter 9 of (Bishop, 2006).

Selecting the number of clusters
---------------------------------------------

The number of clusters is not known beforehand. More mixture components always leads to a higher likelihood score. Therefore, to be able to use split and delete steps, and sensibly determine the number of clusters, we must add a penalty term based on number of clusters--a penalty for model complexity. The penalty we use is based on the Bayes Information Criteria (BIC).

To understand the BIC we must introduce the Bayesian perspective on inference (as opposed to the maximum likelihood approach). In the Bayesian approach, we are looking for the posterior probabilities like P(Model | Data) and P(Parameters | Data). This requires us to assign prior probabilities P(Model) and P(Parameters).

Terminology note: 

   *  "Parameter" refers to numerical values, like the cluster means :math:`\mu`.
   * "Model" refers to the number of mixture components, i.e., we have a 5-cluster model and 6-cluster model for the data.

The Bayesian approach gives a sensible way to do model selection, i.e., determine the number of clusters. We simply look for the model that has the highest posterior probability P(Model | Data).

   P(Model | Data) = P(Data | Model) P(Model) / P(Data)
   
Ignoring P(Data)--which is summed over all models--and assuming the same prior probability P(Model) for all all models, we are left with:

   P(Model | Data) ~ P(Data | Model)
   
P(Data | Model) is really an integral over all the possible parameters that go with the model. Using the notation X = Data, :math:`\Theta` = Parameters, :math:`\mathcal{M}` = Model.

.. math::
   :nowrap:

   \begin{align*}
      P(X | \mathcal{M}) &= \int d\Theta P(X | \Theta,\mathcal{M}) P(\Theta,\mathcal{M}) \\
      &\approx P(X | \Theta_{MAP},\mathcal{M}) V_{\text{posterior}} \frac{1}{V_{\text{prior}}}
   \end{align*}
   
:math:`V_{\text{posterior}}` is the volume of the peak in the posterior distribution over parameters :math:`P(\Theta|X,\mathcal{M})`, and :math:`V_{\text{prior}}` is the volume of the peak in the prior distribution over parameters :math:`P(\Theta|\mathcal{M})`. :math:`\Theta_{MAP}` is the "maximum a posteriori" value of :math:`\Theta`, i.e., the value that maximizes :math:`P(X | \Theta_{MAP})`. Here we used the approximation of a Gaussian integral as the peak value times its width. Thus :math:`P(\Theta,\mathcal{M}) \approx \frac{1}{V_{\text{prior}}}` and :math:`P(X | \Theta,\mathcal{M}) \approx P(X | \Theta_{MAP},\mathcal{M}) V_{\text{posterior}}`.

The term :math:`\frac{V_{\text{posterior}}}{V_{\text{prior}}}` acts like a "complexity penalty" that penalizes models with more clusters. Assuming scalar parameters :math:`w_1,w_2,\dots,w_p`, then 

.. math:: `\prod_{i=1}^p \frac{\delta w_{i,\text{prior}}}{\delta w_{i,\text{posterior}}}=\left(\frac{\delta w_{\text{prior}}}{\delta w_{\text{posterior}}}\right)^p`

Thus we see that this term decreases exponentially with number of parameters. 

The idea behind the BIC is to try to approximate the Bayesian calculation of the posterior probability :math:`P(X | \mathcal{M}) = P(X | \Theta_{MAP},\mathcal{M}) \frac{V_{\text{posterior}}}{V_{\text{prior}}}`. Then we find the model that maximizes this. :math:`P(X | \Theta_{MAP}) \approx P(X | \Theta_{ML})`. So we have to make a rough approximation of :math:`\frac{V_{\text{posterior}}}{V_{\text{prior}}}`, even though we do not have a prior or a posterior distribution.

For a given scalar parameter :math:`w`, suppose the width of the prior is :math:`\delta w_{\text{prior}}`. Then, after N measurements, we get :math:`\delta w_{\text{posterior}} \approx \delta w_{\text{prior}}/\sqrt N`. Using the above equation, we get

.. math:: 
   :nowrap:

   \begin{align*}
	\frac{V_{\text{posterior}}}{V_{\text{prior}}} &\approx \left( \frac{1}{\sqrt N} \right)^p \\
	\log \frac{V_{\text{posterior}}}{V_{\text{prior}}} &\approx p/2 \log N = BIC
   \end{align*} 

Here N is the number of points, and p is the number of parameters.

Using the BIC, we do model selection by maximizing

.. math:: \log P(X | \Theta_{ML}) - BIC = \log P(X | \Theta_{ML}) - p/2 \log N


In the case of clustering, the parameters of a given cluster are only constrained by the points in that cluster--not all points--so the following is more accurate:

.. math:: BIC = \sum_m p/2 \log N_m

where :math:`N_m` is the number of points in the mth cluster. This amounts to a constant correction of :math:`p/2 \log M`, where M is the total number of clusters. Note that smaller clusters are penalized less than larger clusters, though a small cluster and a large cluster are penalized more than two medium-sized clusters. Interestingly, this resembles the minimum message length, which is also frequently used as a complexity penalty for model selection. In clustering, it takes the form (Wallace et al., 1987) via (Shoham et al., 2003):

.. math:: MML = \sum_m p/2 \log N_m/12 + M/2 \log N/12 + M(p+1)/2


.. _KlustaKwik: http://sourceforge.org/klustakwik

* Bishop, M. *Pattern recognition and machine learning* (2006).
* Shoham, S. and Fellows, M.R. and Normann, R.A. *Robust, automatic spike sorting using mixtures of multivariate t-distributions*, Journal of neuroscience methods (2003).
* Wallace, C.S. and Freeman, P.R. *Estimation and inference by compact coding*. Journal of the royal statistical society (1987).

