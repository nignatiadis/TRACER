function fast_minim(
   rcount;
   max_iter = 7,
   λ = 0.01,
   weighted = true,
   threshold = 10^(-4),
   refit = true,
   P1_sparsity = fill(1, (size(rcount, 1), size(rcount, 1))),
   P2_sparsity = fill(1, (size(rcount, 1), size(rcount, 1))),
)

   obj_vals = ones(Float64, max_iter) .* Inf
   if (weighted)
      wts = 1 ./ max.(rcount, 2.0)
      wts = wts ./ sum(wts) .* length(wts)
   else
      wts = fill(1.0, size(rcount))
   end

   nclust = size(rcount, 1)
   P1 = Matrix(1.0I, nclust, nclust)
   P2 = Matrix(1.0I, nclust, nclust)

   rs = rcount ./ sum(rcount, dims = 1)

   t_final = size(rs, 2)
   K_clusters = size(rs, 1)

   jump_days = [3, 4, 6, 7]

   obj_vals = fill(Inf, max_iter)

   n_transitions = 2

   theta_m = (rs[:, jump_days] .+ rs[:, jump_days.-1]) ./ 2 # initialize
   theta_m_l = theta_m[:, 1:2]
   theta_m_r = theta_m[:, 3:4]

   for k = 1:max_iter

      m1 = Model(solver = GurobiSolver(OutputFlag = 0))

      @variable(m1, theta_pred_l[1:K_clusters, 1:4])
      @variable(m1, theta_pred_r[1:K_clusters, 1:3])

      @variable(m1, transition_P1[1:K_clusters, 1:K_clusters] >= 0)
      @variable(m1, transition_P2[1:K_clusters, 1:K_clusters] >= 0)


      for i = 1:K_clusters
         @constraint(m1, sum(transition_P1[:, i]) == 1)
         @constraint(m1, sum(transition_P2[:, i]) == 1)
         for j = 1:K_clusters
            if (P1_sparsity[i, j] == 0)
               @constraint(m1, transition_P1[i, j] == 0)
            end
            if (P2_sparsity[i, j] == 0)
               @constraint(m1, transition_P2[i, j] == 0)
            end
         end

      end

      @constraint(m1, theta_pred_l[:, 2] .== transition_P1 * rs[:, 1]) #day 0->2
      @constraint(m1, theta_pred_l[:, 3:4] .== transition_P1 * theta_m_l) # day 2->6,6->10
      @constraint(m1, theta_pred_r[:, 1] .== transition_P2 * rs[:, 4])# day 10->w2
      @constraint(m1, theta_pred_r[:, 2:3] .== transition_P2 * theta_m_r) #day w2->w6

      @objective(
         m1,
         Min,
         sum(wts[:, 1:4] .* (theta_pred_l .- rs[:, 1:4]) .^ 2) +
         sum(wts[:, 3:4] .* (theta_m_l .- transition_P1 * rs[:, 2:3]) .^ 2) +
         2 * λ / n_transitions * (sum(transition_P1) - tr(transition_P1)) +
         sum(wts[:, 5:7] .* (theta_pred_r .- rs[:, 5:7]) .^ 2) +
         sum(wts[:, 6:7] .* (theta_m_r - transition_P2 * rs[:, 5:6]) .^ 2) +
         2 * λ / n_transitions * (sum(transition_P2) - tr(transition_P2)),
      )

      status = solve(m1)

      obj_vals[k] = getobjectivevalue(m1)

      P1 = getvalue(transition_P1)
      P2 = getvalue(transition_P2)

      m2 = Model(solver = GurobiSolver(OutputFlag = 0))
      @variable(m2, theta_m_l_c[1:K_clusters, 1:2] >= 0)
      @variable(m2, theta_m_r_c[1:K_clusters, 1:2] >= 0)
      @variable(m2, theta_pred_l_c[1:K_clusters, 1:4])
      @variable(m2, theta_pred_r_c[1:K_clusters, 1:3])

      for i = 1:2
         @constraint(m2, sum(theta_m_l_c[:, i]) == 1)
         @constraint(m2, sum(theta_m_r_c[:, i]) == 1)
      end

      @constraint(m2, theta_pred_l_c[:, 2] .== P1 * rs[:, 1]) #day 0->2
      @constraint(m2, theta_pred_l_c[:, 3:4] .== P1 * theta_m_l_c) # day 2->6,6->10
      @constraint(m2, theta_pred_r_c[:, 1] .== P2 * rs[:, 4])# day 10->w2
      @constraint(m2, theta_pred_r_c[:, 2:3] .== P2 * theta_m_r_c) #day w2->w6

      @objective(
         m2,
         Min,
         sum(wts[:, 1:4] .* (theta_pred_l_c .- rs[:, 1:4]) .^ 2) +
         sum(wts[:, 3:4] .* (theta_m_l_c .- P1 * rs[:, 2:3]) .^ 2) +
         sum(wts[:, 5:7] .* (theta_pred_r_c .- rs[:, 5:7]) .^ 2) +
         sum(wts[:, 6:7] .* (theta_m_r_c - P2 * rs[:, 5:6]) .^ 2),
      )

      status2 = solve(m2)

      theta_m_l = getvalue(theta_m_l_c)
      theta_m_r = getvalue(theta_m_r_c)

   end

   P1_sparsity = 1 .* (P1 .> threshold)
   P2_sparsity = 1 .* (P2 .> threshold)

   res = Dict(
      "P1" => P1,
      "P2" => P2,
      "obj_vals" => obj_vals,
      "theta_m" => [theta_m_l theta_m_r],
   )

   if refit
      post_res = fast_minim(
         rcount;
         max_iter = max_iter,
         λ = 0.0,
         weighted = weighted,
         threshold = threshold,
         refit = false,
         P1_sparsity = P1_sparsity,
         P2_sparsity = P2_sparsity,
      )

      res = Dict("res" => res, "post_res" => post_res)
   end

   res
end
