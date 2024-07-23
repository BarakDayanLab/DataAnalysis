import inspect
from qm.qua import *
from Experiments.Enums.IOParameters import IOParameters as IOP


# -------------------------------------------------------------------
# OPX Stuff - can be used from different OPX_Code in all experiments
# -------------------------------------------------------------------


class OPX_Utils:

    def __init__(self):
        raise Exception('No need to initialize OPX_Utils class. It used for statis methods only')

    @staticmethod
    def assign_experiment_variables(self_obj):
        """
        Populate the parameters with default values
        """
        MAX_SIZE_OF_PARAMS_VECTOR_IN_OPX = 100

        # Declare a vector of 100 entries for all possible parameters
        params = declare(int, size=MAX_SIZE_OF_PARAMS_VECTOR_IN_OPX)

        parameters_count = 0
        self_members = inspect.getmembers(self_obj)

        # Iterate over all "self" values
        for member in self_members:
            member_name = member[0]
            if member_name.startswith('_'):
                continue

            member_value = member[1]
            member_value_type = type(member_value)
            if member_value_type != int and member_value_type != float and member_value_type != bool:
                continue

            if IOP.has(member_name):
                index = IOP.value_from_string(member_name)
                assign(params[index], member_value)
                parameters_count += 1

        # Iterate over all experiment values
        for key, value in self_obj.Exp_Values.items():
            if IOP.has(key):
                index = IOP.value_from_string(key)
                assign(params[index], value)
                parameters_count += 1

        return params

    @staticmethod
    def parameters_update(param):
        return

    # TODO: Q:
    # TODO: Need to check what this does. It seems to do nothing as it multiplies by ZERO.
    # TODO: Dor says it's a hack by Quantum Machine that is suppose to tie variables to elements
    # TODO: Need to check it.
    #
    @staticmethod
    def assign_variables_to_element(element, *variables):
        """
        Forces the given variables to be used by the given element thread. Useful as a workaround for when the compiler
        wrongly assigns variables which can cause gaps.
        To be used at the beginning of a program, will add a 16ns wait to the given element. Use an `align()` if needed.

        Example::

            >>> with program() as program_name:
            >>>     a = declare(int)
            >>>     b = declare(fixed)
            >>>     assign_variables_to_element('resonator', a, b)
            >>>     align()
            >>>     ...

        """
        _exp = variables[0]
        for variable in variables[1:]:
            _exp += variable
        wait(4 + 0 * _exp, element)
