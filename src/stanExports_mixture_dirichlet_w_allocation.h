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
namespace model_mixture_dirichlet_w_allocation_namespace {
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
    reader.add_event(0, 0, "start", "model_mixture_dirichlet_w_allocation");
    reader.add_event(90, 88, "end", "model_mixture_dirichlet_w_allocation");
    return reader;
}
#include <stan_meta_header.hpp>
class model_mixture_dirichlet_w_allocation : public prob_grad {
private:
    int n_groups;
    int n_data;
    vector<vector<double> > y;
    double alpha0;
    double period_over_2pi;
public:
    model_mixture_dirichlet_w_allocation(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_mixture_dirichlet_w_allocation(stan::io::var_context& context__,
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
        static const char* function__ = "model_mixture_dirichlet_w_allocation_namespace::model_mixture_dirichlet_w_allocation";
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
            context__.validate_dims("data initialization", "alpha0", "double", context__.to_vec());
            alpha0 = double(0);
            vals_r__ = context__.vals_r("alpha0");
            pos__ = 0;
            alpha0 = vals_r__[pos__++];
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "period_over_2pi", "double", context__.to_vec());
            period_over_2pi = double(0);
            vals_r__ = context__.vals_r("period_over_2pi");
            pos__ = 0;
            period_over_2pi = vals_r__[pos__++];
            // validate, data variables
            current_statement_begin__ = 5;
            current_statement_begin__ = 6;
            current_statement_begin__ = 7;
            current_statement_begin__ = 8;
            check_greater_or_equal(function__,"alpha0",alpha0,0);
            current_statement_begin__ = 9;
            check_greater_or_equal(function__,"period_over_2pi",period_over_2pi,0);
            // initialize data variables
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 13;
            validate_non_negative_index("mu", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 14;
            validate_non_negative_index("sigma", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 16;
            validate_non_negative_index("nu", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 17;
            validate_non_negative_index("kappa", "n_groups", n_groups);
            num_params_r__ += n_groups;
            current_statement_begin__ = 19;
            validate_non_negative_index("v", "n_groups", n_groups);
            num_params_r__ += n_groups;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_mixture_dirichlet_w_allocation() { }
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
        if (!(context__.contains_r("v")))
            throw std::runtime_error("variable v missing");
        vals_r__ = context__.vals_r("v");
        pos__ = 0U;
        validate_non_negative_index("v", "n_groups", n_groups);
        context__.validate_dims("initialization", "v", "double", context__.to_vec(n_groups));
        std::vector<double> v(n_groups,double(0));
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            v[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < n_groups; ++i0__)
            try {
            writer__.scalar_lub_unconstrain(0,1,v[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable v: ") + e.what());
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
            vector<local_scalar_t__> v;
            size_t dim_v_0__ = n_groups;
            v.reserve(dim_v_0__);
            for (size_t k_0__ = 0; k_0__ < dim_v_0__; ++k_0__) {
                if (jacobian__)
                    v.push_back(in__.scalar_lub_constrain(0,1,lp__));
                else
                    v.push_back(in__.scalar_lub_constrain(0,1));
            }
            // transformed parameters
            current_statement_begin__ = 23;
            validate_non_negative_index("pmix", "n_groups", n_groups);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  pmix(static_cast<Eigen::VectorXd::Index>(n_groups));
            (void) pmix;  // dummy to suppress unused var warning
            stan::math::initialize(pmix, DUMMY_VAR__);
            stan::math::fill(pmix,DUMMY_VAR__);
            current_statement_begin__ = 24;
            stan::model::assign(pmix, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        get_base1(v,1,"v",1), 
                        "assigning variable pmix");
            current_statement_begin__ = 26;
            for (int j = 2; j <= (n_groups - 1); ++j) {
                current_statement_begin__ = 27;
                stan::model::assign(pmix, 
                            stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                            stan::model::deep_copy((((get_base1(v,j,"v",1) * (1 - get_base1(v,(j - 1),"v",1))) * get_base1(pmix,(j - 1),"pmix",1)) / get_base1(v,(j - 1),"v",1))), 
                            "assigning variable pmix");
            }
            current_statement_begin__ = 29;
            stan::model::assign(pmix, 
                        stan::model::cons_list(stan::model::index_uni(n_groups), stan::model::nil_index_list()), 
                        stan::model::deep_copy((1 - sum(stan::model::rvalue(pmix, stan::model::cons_list(stan::model::index_min_max(1, (n_groups - 1)), stan::model::nil_index_list()), "pmix")))), 
                        "assigning variable pmix");
            // validate transformed parameters
            for (int i0__ = 0; i0__ < n_groups; ++i0__) {
                if (stan::math::is_uninitialized(pmix(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: pmix" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 23;
            stan::math::check_simplex(function__,"pmix",pmix);
            // model body
            {
            current_statement_begin__ = 33;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            vector<local_scalar_t__> contributions(n_groups);
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions,DUMMY_VAR__);
            current_statement_begin__ = 36;
            lp_accum__.add(cauchy_log<propto__>(mu, 0, 10));
            current_statement_begin__ = 37;
            lp_accum__.add(cauchy_log<propto__>(sigma, 0, 10));
            current_statement_begin__ = 39;
            lp_accum__.add(cauchy_log<propto__>(nu, 0, 10));
            current_statement_begin__ = 40;
            lp_accum__.add(cauchy_log<propto__>(kappa, 0, 10));
            current_statement_begin__ = 42;
            lp_accum__.add(beta_log<propto__>(v, 1, alpha0));
            current_statement_begin__ = 46;
            for (int i = 1; i <= n_data; ++i) {
                current_statement_begin__ = 47;
                for (int k = 1; k <= n_groups; ++k) {
                    current_statement_begin__ = 48;
                    if (as_bool(logical_lt(get_base1(kappa,k,"kappa",1),100))) {
                        current_statement_begin__ = 49;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + von_mises_log((get_base1(get_base1(y,i,"y",1),2,"y",2) / period_over_2pi),get_base1(nu,k,"nu",1),get_base1(kappa,k,"kappa",1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 55;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + normal_log((get_base1(get_base1(y,i,"y",1),2,"y",2) / period_over_2pi),get_base1(nu,k,"nu",1),(1 / stan::math::sqrt(get_base1(kappa,k,"kappa",1))))), 
                                    "assigning variable contributions");
                    }
                }
                current_statement_begin__ = 60;
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
        names__.push_back("v");
        names__.push_back("pmix");
        names__.push_back("log_lik");
        names__.push_back("contributions");
        names__.push_back("post_prob");
        names__.push_back("allocation");
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
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_data);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_data);
        dims__.push_back(n_groups);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_data);
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
        static const char* function__ = "model_mixture_dirichlet_w_allocation_namespace::write_array";
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
        vector<double> v;
        size_t dim_v_0__ = n_groups;
        for (size_t k_0__ = 0; k_0__ < dim_v_0__; ++k_0__) {
            v.push_back(in__.scalar_lub_constrain(0,1));
        }
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
            vars__.push_back(v[k_0__]);
            }
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            current_statement_begin__ = 23;
            validate_non_negative_index("pmix", "n_groups", n_groups);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  pmix(static_cast<Eigen::VectorXd::Index>(n_groups));
            (void) pmix;  // dummy to suppress unused var warning
            stan::math::initialize(pmix, DUMMY_VAR__);
            stan::math::fill(pmix,DUMMY_VAR__);
            current_statement_begin__ = 24;
            stan::model::assign(pmix, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        get_base1(v,1,"v",1), 
                        "assigning variable pmix");
            current_statement_begin__ = 26;
            for (int j = 2; j <= (n_groups - 1); ++j) {
                current_statement_begin__ = 27;
                stan::model::assign(pmix, 
                            stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                            stan::model::deep_copy((((get_base1(v,j,"v",1) * (1 - get_base1(v,(j - 1),"v",1))) * get_base1(pmix,(j - 1),"pmix",1)) / get_base1(v,(j - 1),"v",1))), 
                            "assigning variable pmix");
            }
            current_statement_begin__ = 29;
            stan::model::assign(pmix, 
                        stan::model::cons_list(stan::model::index_uni(n_groups), stan::model::nil_index_list()), 
                        stan::model::deep_copy((1 - sum(stan::model::rvalue(pmix, stan::model::cons_list(stan::model::index_min_max(1, (n_groups - 1)), stan::model::nil_index_list()), "pmix")))), 
                        "assigning variable pmix");
            // validate transformed parameters
            current_statement_begin__ = 23;
            stan::math::check_simplex(function__,"pmix",pmix);
            // write transformed parameters
            if (include_tparams__) {
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(pmix[k_0__]);
            }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 67;
            validate_non_negative_index("log_lik", "n_data", n_data);
            vector<local_scalar_t__> log_lik(n_data);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik,DUMMY_VAR__);
            current_statement_begin__ = 68;
            validate_non_negative_index("contributions", "n_groups", n_groups);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  contributions(static_cast<Eigen::VectorXd::Index>(n_groups));
            (void) contributions;  // dummy to suppress unused var warning
            stan::math::initialize(contributions, DUMMY_VAR__);
            stan::math::fill(contributions,DUMMY_VAR__);
            current_statement_begin__ = 69;
            validate_non_negative_index("post_prob", "n_groups", n_groups);
            validate_non_negative_index("post_prob", "n_data", n_data);
            vector<Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> > post_prob(n_data, (Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> (static_cast<Eigen::VectorXd::Index>(n_groups))));
            stan::math::initialize(post_prob, DUMMY_VAR__);
            stan::math::fill(post_prob,DUMMY_VAR__);
            current_statement_begin__ = 70;
            validate_non_negative_index("allocation", "n_data", n_data);
            vector<int> allocation(n_data, 0);
            stan::math::fill(allocation, std::numeric_limits<int>::min());
            current_statement_begin__ = 72;
            for (int i = 1; i <= n_data; ++i) {
                current_statement_begin__ = 73;
                for (int k = 1; k <= n_groups; ++k) {
                    current_statement_begin__ = 74;
                    if (as_bool(logical_lt(get_base1(kappa,k,"kappa",1),100))) {
                        current_statement_begin__ = 75;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + von_mises_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),get_base1(kappa,k,"kappa",1))), 
                                    "assigning variable contributions");
                    } else {
                        current_statement_begin__ = 79;
                        stan::model::assign(contributions, 
                                    stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                    ((stan::math::log(get_base1(pmix,k,"pmix",1)) + normal_log(get_base1(get_base1(y,i,"y",1),1,"y",2),get_base1(mu,k,"mu",1),get_base1(sigma,k,"sigma",1))) + normal_log(get_base1(get_base1(y,i,"y",1),2,"y",2),get_base1(nu,k,"nu",1),(1 / stan::math::sqrt(get_base1(kappa,k,"kappa",1))))), 
                                    "assigning variable contributions");
                    }
                }
                current_statement_begin__ = 84;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            log_sum_exp(contributions), 
                            "assigning variable log_lik");
                current_statement_begin__ = 85;
                stan::model::assign(post_prob, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            softmax(contributions), 
                            "assigning variable post_prob");
                current_statement_begin__ = 86;
                stan::model::assign(allocation, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            categorical_rng(get_base1(post_prob,i,"post_prob",1), base_rng__), 
                            "assigning variable allocation");
            }
            // validate generated quantities
            current_statement_begin__ = 67;
            current_statement_begin__ = 68;
            current_statement_begin__ = 69;
            for (int k0__ = 0; k0__ < n_data; ++k0__) {
                stan::math::check_simplex(function__,"post_prob[k0__]",post_prob[k0__]);
            }
            current_statement_begin__ = 70;
            // write generated quantities
            for (int k_0__ = 0; k_0__ < n_data; ++k_0__) {
            vars__.push_back(log_lik[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_groups; ++k_0__) {
            vars__.push_back(contributions[k_0__]);
            }
            for (int k_1__ = 0; k_1__ < n_groups; ++k_1__) {
                for (int k_0__ = 0; k_0__ < n_data; ++k_0__) {
                vars__.push_back(post_prob[k_0__][k_1__]);
                }
            }
            for (int k_0__ = 0; k_0__ < n_data; ++k_0__) {
            vars__.push_back(allocation[k_0__]);
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
        return "model_mixture_dirichlet_w_allocation";
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
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pmix" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
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
        for (int k_1__ = 1; k_1__ <= n_groups; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "post_prob" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "allocation" << '.' << k_0__;
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
        for (int k_0__ = 1; k_0__ <= n_groups; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= (n_groups - 1); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pmix" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
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
        for (int k_1__ = 1; k_1__ <= (n_groups - 1); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "post_prob" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_0__ = 1; k_0__ <= n_data; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "allocation" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}
typedef model_mixture_dirichlet_w_allocation_namespace::model_mixture_dirichlet_w_allocation stan_model;
#endif
