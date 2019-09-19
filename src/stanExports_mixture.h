// Generated by rstantools.  Do not edit by hand.

/*
    bayesecopop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bayesecopop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bayesecopop.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_mixture_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_mixture");
    reader.add_event(78, 76, "end", "model_mixture");
    return reader;
}
#include <stan_meta_header.hpp>
class model_mixture : public prob_grad {
private:
        int n_groups;
        int n_data;
        std::vector<std::vector<double> > y;
        vector_d alpha;
public:
    model_mixture(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_mixture(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_mixture_namespace::model_mixture";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "n_groups", "int", context__.to_vec());
            n_groups = int(0);
            vals_i__ = context__.vals_i("n_groups");
            pos__ = 0;
            n_groups = vals_i__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "n_data", "int", context__.to_vec());
            n_data = int(0);
            vals_i__ = context__.vals_i("n_data");
            pos__ = 0;
            n_data = vals_i__[pos__++];
            current_statement_begin__ = 7;
            validate_non_negative_index("y", "n_data", n_data);
            validate_non_negative_index("y", "2", 2);
            context__.validate_dims("data initialization", "y", "double", context__.to_vec(n_data,2));
            y = std::vector<std::vector<double> >(n_data, std::vector<double>(2, double(0)));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_k_0_max__ = n_data;
            size_t y_k_1_max__ = 2;
            for (size_t k_1__ = 0; k_1__ < y_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                    y[k_0__][k_1__] = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("alpha", "n_groups", n_groups);
            context__.validate_dims("data initialization", "alpha", "vector_d", context__.to_vec(n_groups));
            alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>(n_groups);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            size_t alpha_j_1_max__ = n_groups;
            for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
                alpha(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 12;
            validate_non_negative_index("mu", "n_groups", n_groups);
            num_params_r__ += (1 * n_groups);
            current_statement_begin__ = 13;
            validate_non_negative_index("sigma", "n_groups", n_groups);
            num_params_r__ += (1 * n_groups);
            current_statement_begin__ = 15;
            validate_non_negative_index("nu", "n_groups", n_groups);
            num_params_r__ += (1 * n_groups);
            current_statement_begin__ = 16;
            validate_non_negative_index("kappa", "n_groups", n_groups);
            num_params_r__ += (1 * n_groups);
            current_statement_begin__ = 18;
            validate_non_negative_index("pmix", "n_groups", n_groups);
            num_params_r__ += (n_groups - 1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_mixture() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 12;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "n_groups", n_groups);
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec(n_groups));
        std::vector<double> mu(n_groups, double(0));
        size_t mu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
            mu[k_0__] = vals_r__[pos__++];
        }
        size_t mu_i_0_max__ = n_groups;
        for (size_t i_0__ = 0; i_0__ < mu_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(mu[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 13;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        validate_non_negative_index("sigma", "n_groups", n_groups);
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec(n_groups));
        std::vector<double> sigma(n_groups, double(0));
        size_t sigma_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < sigma_k_0_max__; ++k_0__) {
            sigma[k_0__] = vals_r__[pos__++];
        }
        size_t sigma_i_0_max__ = n_groups;
        for (size_t i_0__ = 0; i_0__ < sigma_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_lb_unconstrain(0, sigma[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 15;
        if (!(context__.contains_r("nu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable nu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("nu");
        pos__ = 0U;
        validate_non_negative_index("nu", "n_groups", n_groups);
        context__.validate_dims("parameter initialization", "nu", "double", context__.to_vec(n_groups));
        std::vector<double> nu(n_groups, double(0));
        size_t nu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < nu_k_0_max__; ++k_0__) {
            nu[k_0__] = vals_r__[pos__++];
        }
        size_t nu_i_0_max__ = n_groups;
        for (size_t i_0__ = 0; i_0__ < nu_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(nu[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable nu: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 16;
        if (!(context__.contains_r("kappa")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable kappa missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("kappa");
        pos__ = 0U;
        validate_non_negative_index("kappa", "n_groups", n_groups);
        context__.validate_dims("parameter initialization", "kappa", "double", context__.to_vec(n_groups));
        std::vector<double> kappa(n_groups, double(0));
        size_t kappa_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < kappa_k_0_max__; ++k_0__) {
            kappa[k_0__] = vals_r__[pos__++];
        }
        size_t kappa_i_0_max__ = n_groups;
        for (size_t i_0__ = 0; i_0__ < kappa_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_lb_unconstrain(0, kappa[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable kappa: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 18;
        if (!(context__.contains_r("pmix")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable pmix missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("pmix");
        pos__ = 0U;
        validate_non_negative_index("pmix", "n_groups", n_groups);
        context__.validate_dims("parameter initialization", "pmix", "vector_d", context__.to_vec(n_groups));
        Eigen::Matrix<double, Eigen::Dynamic, 1> pmix(n_groups);
        size_t pmix_j_1_max__ = n_groups;
        for (size_t j_1__ = 0; j_1__ < pmix_j_1_max__; ++j_1__) {
            pmix(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.simplex_unconstrain(pmix);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable pmix: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 12;
            std::vector<local_scalar_t__> mu;
            size_t mu_d_0_max__ = n_groups;
            mu.reserve(mu_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < mu_d_0_max__; ++d_0__) {
                if (jacobian__)
                    mu.push_back(in__.scalar_constrain(lp__));
                else
                    mu.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 13;
            std::vector<local_scalar_t__> sigma;
            size_t sigma_d_0_max__ = n_groups;
            sigma.reserve(sigma_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < sigma_d_0_max__; ++d_0__) {
                if (jacobian__)
                    sigma.push_back(in__.scalar_lb_constrain(0, lp__));
                else
                    sigma.push_back(in__.scalar_lb_constrain(0));
            }
            current_statement_begin__ = 15;
            std::vector<local_scalar_t__> nu;
            size_t nu_d_0_max__ = n_groups;
            nu.reserve(nu_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < nu_d_0_max__; ++d_0__) {
                if (jacobian__)
                    nu.push_back(in__.scalar_constrain(lp__));
                else
                    nu.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 16;
            std::vector<local_scalar_t__> kappa;
            size_t kappa_d_0_max__ = n_groups;
            kappa.reserve(kappa_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < kappa_d_0_max__; ++d_0__) {
                if (jacobian__)
                    kappa.push_back(in__.scalar_lb_constrain(0, lp__));
                else
                    kappa.push_back(in__.scalar_lb_constrain(0));
            }
            current_statement_begin__ = 18;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> pmix;
            (void) pmix;  // dummy to suppress unused var warning
            if (jacobian__)
                pmix = in__.simplex_constrain(n_groups, lp__);
            else
                pmix = in__.simplex_constrain(n_groups);
            // model body
            {
            current_statement_begin__ = 22;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            std::vector<local_scalar_t__  > contributions(n_groups, local_scalar_t__(DUMMY_VAR__));
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions, DUMMY_VAR__);
            current_statement_begin__ = 25;
            lp_accum__.add(cauchy_log<propto__>(mu, 0, 10));
            current_statement_begin__ = 26;
            lp_accum__.add(cauchy_log<propto__>(sigma, 0, 10));
            current_statement_begin__ = 28;
            lp_accum__.add(cauchy_log<propto__>(nu, 0, 10));
            current_statement_begin__ = 29;
            lp_accum__.add(cauchy_log<propto__>(kappa, 0, 10));
            current_statement_begin__ = 31;
            lp_accum__.add(dirichlet_log<propto__>(pmix, alpha));
            current_statement_begin__ = 35;
            for (int i = 1; i <= n_data; ++i) {
                current_statement_begin__ = 36;
                for (int k = 1; k <= n_groups; ++k) {
                    current_statement_begin__ = 37;
                    if (as_bool(logical_lt(get_base1(kappa, k, "kappa", 1), 100))) {
                        current_statement_begin__ = 38;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix, k, "pmix", 1)) + normal_log(get_base1(get_base1(y, i, "y", 1), 1, "y", 2), get_base1(mu, k, "mu", 1), get_base1(sigma, k, "sigma", 1))) + von_mises_log(get_base1(get_base1(y, i, "y", 1), 2, "y", 2), get_base1(nu, k, "nu", 1), get_base1(kappa, k, "kappa", 1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 44;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix, k, "pmix", 1)) + normal_log(get_base1(get_base1(y, i, "y", 1), 1, "y", 2), get_base1(mu, k, "mu", 1), get_base1(sigma, k, "sigma", 1))) + normal_log(get_base1(get_base1(y, i, "y", 1), 2, "y", 2), get_base1(nu, k, "nu", 1), (1 / stan::math::sqrt(get_base1(kappa, k, "kappa", 1))))), 
                                    "assigning variable contributions");
                    }
                }
                current_statement_begin__ = 49;
                lp_accum__.add(log_sum_exp(contributions));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu");
        names__.push_back("sigma");
        names__.push_back("nu");
        names__.push_back("kappa");
        names__.push_back("pmix");
        names__.push_back("log_lik");
        names__.push_back("contributions");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_data);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_mixture_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        std::vector<double> mu;
        size_t mu_d_0_max__ = n_groups;
        mu.reserve(mu_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < mu_d_0_max__; ++d_0__) {
            mu.push_back(in__.scalar_constrain());
        }
        size_t mu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
            vars__.push_back(mu[k_0__]);
        }
        std::vector<double> sigma;
        size_t sigma_d_0_max__ = n_groups;
        sigma.reserve(sigma_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < sigma_d_0_max__; ++d_0__) {
            sigma.push_back(in__.scalar_lb_constrain(0));
        }
        size_t sigma_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < sigma_k_0_max__; ++k_0__) {
            vars__.push_back(sigma[k_0__]);
        }
        std::vector<double> nu;
        size_t nu_d_0_max__ = n_groups;
        nu.reserve(nu_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < nu_d_0_max__; ++d_0__) {
            nu.push_back(in__.scalar_constrain());
        }
        size_t nu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < nu_k_0_max__; ++k_0__) {
            vars__.push_back(nu[k_0__]);
        }
        std::vector<double> kappa;
        size_t kappa_d_0_max__ = n_groups;
        kappa.reserve(kappa_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < kappa_d_0_max__; ++d_0__) {
            kappa.push_back(in__.scalar_lb_constrain(0));
        }
        size_t kappa_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < kappa_k_0_max__; ++k_0__) {
            vars__.push_back(kappa[k_0__]);
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> pmix = in__.simplex_constrain(n_groups);
        size_t pmix_j_1_max__ = n_groups;
        for (size_t j_1__ = 0; j_1__ < pmix_j_1_max__; ++j_1__) {
            vars__.push_back(pmix(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 55;
            validate_non_negative_index("log_lik", "n_data", n_data);
            std::vector<double> log_lik(n_data, double(0));
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            current_statement_begin__ = 56;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            std::vector<double> contributions(n_groups, double(0));
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 60;
            for (int i = 1; i <= n_data; ++i) {
                current_statement_begin__ = 61;
                for (int k = 1; k <= n_groups; ++k) {
                    current_statement_begin__ = 62;
                    if (as_bool(logical_lt(get_base1(kappa, k, "kappa", 1), 100))) {
                        current_statement_begin__ = 63;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix, k, "pmix", 1)) + normal_log(get_base1(get_base1(y, i, "y", 1), 1, "y", 2), get_base1(mu, k, "mu", 1), get_base1(sigma, k, "sigma", 1))) + von_mises_log(get_base1(get_base1(y, i, "y", 1), 2, "y", 2), get_base1(nu, k, "nu", 1), get_base1(kappa, k, "kappa", 1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 67;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix, k, "pmix", 1)) + normal_log(get_base1(get_base1(y, i, "y", 1), 1, "y", 2), get_base1(mu, k, "mu", 1), get_base1(sigma, k, "sigma", 1))) + normal_log(get_base1(get_base1(y, i, "y", 1), 2, "y", 2), get_base1(nu, k, "nu", 1), (1 / stan::math::sqrt(get_base1(kappa, k, "kappa", 1))))), 
                                    "assigning variable contributions");
                    }
                }
                current_statement_begin__ = 72;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            log_sum_exp(contributions), 
                            "assigning variable log_lik");
            }
            // validate, write generated quantities
            current_statement_begin__ = 55;
            size_t log_lik_k_0_max__ = n_data;
            for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
                vars__.push_back(log_lik[k_0__]);
            }
            current_statement_begin__ = 56;
            size_t contributions_k_0_max__ = n_groups;
            for (size_t k_0__ = 0; k_0__ < contributions_k_0_max__; ++k_0__) {
                vars__.push_back(contributions[k_0__]);
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_mixture";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sigma_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < sigma_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t nu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < nu_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nu" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t kappa_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < kappa_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "kappa" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t pmix_j_1_max__ = n_groups;
        for (size_t j_1__ = 0; j_1__ < pmix_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pmix" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t log_lik_k_0_max__ = n_data;
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t contributions_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < contributions_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "contributions" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < mu_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sigma_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < sigma_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t nu_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < nu_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nu" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t kappa_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < kappa_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "kappa" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t pmix_j_1_max__ = (n_groups - 1);
        for (size_t j_1__ = 0; j_1__ < pmix_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pmix" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t log_lik_k_0_max__ = n_data;
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t contributions_k_0_max__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < contributions_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "contributions" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_mixture_namespace::model_mixture stan_model;
#endif