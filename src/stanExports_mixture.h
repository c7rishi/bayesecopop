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
// Code generated by Stan version 2.18.1
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
    vector<vector<double> > y;
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
        // initialize member variables
        try {
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
            validate_non_negative_index("y", "n_data", n_data);
            validate_non_negative_index("y", "2", 2);
            y = std::vector<std::vector<double> >(n_data,std::vector<double>(2,double(0)));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_limit_1__ = 2;
            for (size_t i_1__ = 0; i_1__ < y_limit_1__; ++i_1__) {
                size_t y_limit_0__ = n_data;
                for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
                    y[i_0__][i_1__] = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("alpha", "n_groups", n_groups);
            context__.validate_dims("data initialization", "alpha", "vector_d", context__.to_vec(n_groups));
            validate_non_negative_index("alpha", "n_groups", n_groups);
            alpha = vector_d(static_cast<Eigen::VectorXd::Index>(n_groups));
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            size_t alpha_i_vec_lim__ = n_groups;
            for (size_t i_vec__ = 0; i_vec__ < alpha_i_vec_lim__; ++i_vec__) {
                alpha[i_vec__] = vals_r__[pos__++];
            }
            // validate, data variables
            current_statement_begin__ = 5;
            current_statement_begin__ = 6;
            current_statement_begin__ = 7;
            current_statement_begin__ = 8;
            // initialize data variables
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 12;
            validate_non_negative_index("mu", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 13;
            validate_non_negative_index("sigma", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 15;
            validate_non_negative_index("nu", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 16;
            validate_non_negative_index("kappa", "n_groups", n_groups);
            num_params_r__ += n_groups;
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
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        if (!(context__.contains_r("mu")))
            throw std::runtime_error("variable mu missing");
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "n_groups", n_groups);
        context__.validate_dims("initialization", "mu", "double", context__.to_vec(n_groups));
        std::vector<double> mu(n_groups,double(0));
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            mu[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            try {
            writer__.scalar_unconstrain(mu[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable mu: ") + e.what());
        }
        if (!(context__.contains_r("sigma")))
            throw std::runtime_error("variable sigma missing");
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        validate_non_negative_index("sigma", "n_groups", n_groups);
        context__.validate_dims("initialization", "sigma", "double", context__.to_vec(n_groups));
        std::vector<double> sigma(n_groups,double(0));
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            sigma[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            try {
            writer__.scalar_lb_unconstrain(0,sigma[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma: ") + e.what());
        }
        if (!(context__.contains_r("nu")))
            throw std::runtime_error("variable nu missing");
        vals_r__ = context__.vals_r("nu");
        pos__ = 0U;
        validate_non_negative_index("nu", "n_groups", n_groups);
        context__.validate_dims("initialization", "nu", "double", context__.to_vec(n_groups));
        std::vector<double> nu(n_groups,double(0));
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            nu[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            try {
            writer__.scalar_unconstrain(nu[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable nu: ") + e.what());
        }
        if (!(context__.contains_r("kappa")))
            throw std::runtime_error("variable kappa missing");
        vals_r__ = context__.vals_r("kappa");
        pos__ = 0U;
        validate_non_negative_index("kappa", "n_groups", n_groups);
        context__.validate_dims("initialization", "kappa", "double", context__.to_vec(n_groups));
        std::vector<double> kappa(n_groups,double(0));
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            kappa[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            try {
            writer__.scalar_lb_unconstrain(0,kappa[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable kappa: ") + e.what());
        }
        if (!(context__.contains_r("pmix")))
            throw std::runtime_error("variable pmix missing");
        vals_r__ = context__.vals_r("pmix");
        pos__ = 0U;
        validate_non_negative_index("pmix", "n_groups", n_groups);
        context__.validate_dims("initialization", "pmix", "vector_d", context__.to_vec(n_groups));
        vector_d pmix(static_cast<Eigen::VectorXd::Index>(n_groups));
        for (int j1__ = 0U; j1__ < n_groups; ++j1__)
            pmix(j1__) = vals_r__[pos__++];
        try {
            writer__.simplex_unconstrain(pmix);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable pmix: ") + e.what());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
            vector<local_scalar_t__> mu;
            size_t dim_mu_0__ = n_groups;
            mu.reserve(dim_mu_0__);
            for (size_t k_0__ = 0; k_0__ < dim_mu_0__; ++k_0__) {
                if (jacobian__)
                    mu.push_back(in__.scalar_constrain(lp__));
                else
                    mu.push_back(in__.scalar_constrain());
            }
            vector<local_scalar_t__> sigma;
            size_t dim_sigma_0__ = n_groups;
            sigma.reserve(dim_sigma_0__);
            for (size_t k_0__ = 0; k_0__ < dim_sigma_0__; ++k_0__) {
                if (jacobian__)
                    sigma.push_back(in__.scalar_lb_constrain(0,lp__));
                else
                    sigma.push_back(in__.scalar_lb_constrain(0));
            }
            vector<local_scalar_t__> nu;
            size_t dim_nu_0__ = n_groups;
            nu.reserve(dim_nu_0__);
            for (size_t k_0__ = 0; k_0__ < dim_nu_0__; ++k_0__) {
                if (jacobian__)
                    nu.push_back(in__.scalar_constrain(lp__));
                else
                    nu.push_back(in__.scalar_constrain());
            }
            vector<local_scalar_t__> kappa;
            size_t dim_kappa_0__ = n_groups;
            kappa.reserve(dim_kappa_0__);
            for (size_t k_0__ = 0; k_0__ < dim_kappa_0__; ++k_0__) {
                if (jacobian__)
                    kappa.push_back(in__.scalar_lb_constrain(0,lp__));
                else
                    kappa.push_back(in__.scalar_lb_constrain(0));
            }
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  pmix;
            (void) pmix;  // dummy to suppress unused var warning
            if (jacobian__)
                pmix = in__.simplex_constrain(n_groups,lp__);
            else
                pmix = in__.simplex_constrain(n_groups);
            // transformed parameters
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // model body
            {
            current_statement_begin__ = 22;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            vector<local_scalar_t__> contributions(n_groups);
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions,DUMMY_VAR__);
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
                    if (as_bool(logical_lt(get_base1(kappa,k,"kappa",1),100))) {
                        current_statement_begin__ = 38;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + von_mises_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),get_base1(kappa,k,"kappa",1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 44;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + normal_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),(1 / stan::math::sqrt(get_base1(kappa,k,"kappa",1))))), 
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
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_mixture_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector<double> mu;
        size_t dim_mu_0__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < dim_mu_0__; ++k_0__) {
            mu.push_back(in__.scalar_constrain());
        }
        vector<double> sigma;
        size_t dim_sigma_0__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < dim_sigma_0__; ++k_0__) {
            sigma.push_back(in__.scalar_lb_constrain(0));
        }
        vector<double> nu;
        size_t dim_nu_0__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < dim_nu_0__; ++k_0__) {
            nu.push_back(in__.scalar_constrain());
        }
        vector<double> kappa;
        size_t dim_kappa_0__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < dim_kappa_0__; ++k_0__) {
            kappa.push_back(in__.scalar_lb_constrain(0));
        }
        vector_d pmix = in__.simplex_constrain(n_groups);
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(mu[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(sigma[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(nu[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(kappa[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(pmix[k_0__]);
            }
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // validate transformed parameters
            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 55;
            validate_non_negative_index("log_lik", "n_data", n_data);
            vector<local_scalar_t__> log_lik(n_data);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik,DUMMY_VAR__);
            current_statement_begin__ = 56;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            vector<local_scalar_t__> contributions(n_groups);
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions,DUMMY_VAR__);
            current_statement_begin__ = 60;
            for (int i = 1; i <= n_data; ++i) {
                current_statement_begin__ = 61;
                for (int k = 1; k <= n_groups; ++k) {
                    current_statement_begin__ = 62;
                    if (as_bool(logical_lt(get_base1(kappa,k,"kappa",1),100))) {
                        current_statement_begin__ = 63;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + von_mises_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),get_base1(kappa,k,"kappa",1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 67;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + normal_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),(1 / stan::math::sqrt(get_base1(kappa,k,"kappa",1))))), 
                                    "assigning variable contributions");
                    }
                }
                current_statement_begin__ = 72;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            log_sum_exp(contributions), 
                            "assigning variable log_lik");
            }
            // validate generated quantities
            current_statement_begin__ = 55;
            current_statement_begin__ = 56;
            // write generated quantities
            for (int k_0__ = 0; k_0__ < n_data; ++k_0__) {
            vars__.push_back(log_lik[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
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
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
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
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "kappa" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pmix" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "contributions" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "kappa" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (n_groups - 1); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pmix" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "contributions" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}
typedef model_mixture_namespace::model_mixture stan_model;
#endif
