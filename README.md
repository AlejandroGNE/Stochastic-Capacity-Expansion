# Stochastic-Capacity-Expansion
 Capacity expansion model with built-in accounting of uncertainty.

This project is under development. Stay tuned!

Electricity system capacity expansion models tend to rely on deterministic inputs, assuming perfect information. Uncertainty is treated through scenario analysis, by re-running the model with alternative sets of inputs. What if we account for uncertainty as part of the optimization, letting it affect the outcome?

This stochastic capacity expansion model achieves that by considering stochastic rather than deterministic inputs.

# Run the example case

This model is written in MATLAB. You will need the following packages: (TODO)

A ficticious example is included. Start from the script 'RunScenario.m', where you can specify a model configuration.

TODO:
- clarify packages/add-ons/toolboxes needed
- explain available model configurations
- switch from mat to csv files so they are visible
- how to run a custom case
    - create new 'RunScenario.m' file with custom name
    - change configurations
    - update model inputs
        - load data
        - existing generators
        - fuel prices
- how to visualize outputs
- host a matlab online server to test?
