#include "OnlinePolicyGradient.h"

PolicyGradientActorCritic::PolicyGradientActorCritic(int n_states_, int n_actions_, real gamma_, real step_size_) :
	n_states(n_states_),
	n_actions(n_actions_),
	gamma(gamma_),
	step_size(step_size_),
	critic(n_states, n_actions, gamma, 0.0, step_size),
	policy(n_states, n_actions),
	params(n_states, n_actions),
	Q(n_states, n_actions)
{
	Reset();
}

real PolicyGradientActorCritic::Observe(real reward, int next_state, int next_action)
{
	if (state >= 0) {
		real d = GradientUpdate(state, action); // not sure how much difference it makes to update the next state-action pair instead
		//printf("%d %d %d %d%f\n", state, action, d);
		Q(state, action) += step_size * (reward + gamma * Q(next_state, next_action) - Q(state, action));
	}
	critic.Observe(reward, next_state, next_action);
	UpdatePolicy();
	state = next_state;
	action = next_action;
	return 0;
}

int PolicyGradientActorCritic::Act(real reward, int next_state)
{
	int next_action = policy.SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}

real PolicyGradientActorCritic::GradientUpdate(int s, int a)
{
	//real U = critic.getValue(s);
	//real U = critic.getValue(s);
	real U = Q(s,a);
	printf("s:%d, a:%d: Q_w:%f Q_S:%f\n", s, a, U, critic.getValue(s,a));
	real delta = 0;
	for (int j=0; j<n_actions; ++j) {
		real p_a = policy.getActionProbability(s, j);
		real d_sj =  0;
		if (j==a) {
			d_sj = (1.0 - p_a) * U;
		} else {
			d_sj = -p_a * U;
		}
		params(s,j) += step_size * d_sj;
		delta += fabs(d_sj);
	}
	return delta;
}

void PolicyGradientActorCritic::UpdatePolicy()
{
	for (int s=0; s<n_states; ++s) {
		Vector eW = exp(params.getRow(s));
		eW /= eW.Sum();
		Vector* pS = policy.getActionProbabilitiesPtr(s);
		(*pS) = eW;
	}
}
