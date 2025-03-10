# Matlab-BHT
 Matlab function for bootstrap-based randomization tests
 
 <div>
  <h1>Test Statistics and Variable Definitions</h1>
  
  <h2>One-sample Test</h2>
  <p>
    <code>T = \frac{\left|\bar{x} - \mu_0\right|}{s/\sqrt{n}}</code>
  </p>
  <ul>
    <li><code>\bar{x}</code>: sample mean</li>
    <li><code>\mu_0</code>: hypothesized (comparison) value</li>
    <li><code>s</code>: sample standard deviation</li>
    <li><code>n</code>: sample size</li>
  </ul>
  
  <h2>Two-sample Test</h2>
  <p>
    <code>T = \frac{\left|\bar{x}_1 - \bar{x}_2\right|}{\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}</code>
  </p>
  <ul>
    <li><code>\bar{x}_1</code> and <code>\bar{x}_2</code>: means of the first and second samples, respectively</li>
    <li><code>s_1^2</code> and <code>s_2^2</code>: sample variances</li>
    <li><code>n_1</code> and <code>n_2</code>: sample sizes for the first and second samples</li>
  </ul>
  
  <h2>Two-sample Paired Test</h2>
  <p>
    <code>T = \frac{\left|\bar{d}\right|}{s_d/\sqrt{n}}</code>
  </p>
  <ul>
    <li><code>\bar{d}</code>: mean of the differences between paired observations</li>
    <li><code>s_d</code>: standard deviation of the differences</li>
    <li><code>n</code>: number of paired observations</li>
  </ul>
  
  <h2>Ranked-consistency Test</h2>
  <p>
    <code>T = \frac{12n}{k(k+1)} \sum_{i=1}^{k} \left(r_i - \frac{k+1}{2}\right)^2</code>
  </p>
  <ul>
    <li><code>n</code>: number of observations (e.g., rows in the data matrix)</li>
    <li><code>k</code>: number of variables (e.g., columns in the data matrix)</li>
    <li><code>r_i</code>: average rank of the <em>i</em>-th variable computed using a tied ranking method</li>
    <li><code>\frac{k+1}{2}</code>: expected rank (midpoint) under the null hypothesis</li>
  </ul>
  
  <h2>Kolmogorov–Smirnov Test</h2>
  <p>
    <code>T = \max_{x} \left| F_1(x) - F_2(x) \right|</code>
  </p>
  <ul>
    <li><code>F_1(x)</code>: empirical cumulative distribution function (ECDF) of the first sample</li>
    <li><code>F_2(x)</code>: ECDF of the second sample</li>
    <li><code>T</code>: maximum absolute difference between the two ECDFs</li>
  </ul>
  
  <h2>Ranksum Test</h2>
  <p>
    <code>T = \left|\sum_{i=1}^{n_1} R_i - \frac{n_1 (n_1+n_2+1)}{2}\right|</code>
  </p>
  <ul>
    <li><code>R_i</code>: rank of the <em>i</em>-th observation from the first sample after combining and ranking both samples</li>
    <li><code>n_1</code>: sample size of the first sample</li>
    <li><code>n_2</code>: sample size of the second sample</li>
    <li><code>\frac{n_1 (n_1+n_2+1)}{2}</code>: expected sum of the ranks for the first sample under the null hypothesis</li>
  </ul>
</div>

