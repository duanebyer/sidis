# Cog script for polarized struct code generation.

import cog
import re

# Which types will be passed by value (instead of by reference).
POD_TYPES = [
    "Real", "float", "double", "long double",
    "bool", "char", "unsigned char",
    "short", "unsigned short", "int", "unsigned",
    "long", "unsigned long", "long long", "unsigned long long",
    "math::Vec3", "math::Bound"]

def mask_find_all(mask, beam_pols, target_pols):
    result = []
    for beam_pol in beam_pols:
        for target_pol in target_pols:
            if mask_match(mask, [beam_pol, target_pol]):
                result.append([beam_pol, target_pol])
    return result

def mask_match(mask, pol):
    match_beam = (mask[0] == "X" or pol[0] == "X" or mask[0] == pol[0])
    match_target = (mask[1] == "X" or pol[1] == "X"
        or (mask[1] == "P" and pol[1] != "U")
        or (pol[1] == "P" and mask[1] != "U")
        or mask[1] == pol[1])
    return match_beam and match_target

def camel_to_snake(name):
    # Taken from `https://stackoverflow.com/questions/1175208`.
    name = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", name).lower()

def indent(string, n=1):
    tabs = "\t" * n
    end = ""
    if len(string) == 0:
        return tabs
    elif string[-1] == "\n":
        end = "\n"
        string = string[:-1]
    result = tabs + string.replace("\n", "\n" + tabs) + end
    return result

def struct_name(base_name, pol):
    return base_name + pol[0] + pol[1]

class Field(object):
    def __init__(self, pol, typ, name):
        self.pol = list(pol)
        self.typ = str(typ)
        self.name = str(name)
    def declare_as_field(self):
        return "{} {}".format(self.typ, self.name)
    def declare_as_arg(self):
        if self.typ in POD_TYPES:
            return "{} {}".format(self.typ, self.name)
        else:
            return "{} const& {}".format(self.typ, self.name)
    def use(self):
        return self.name

class Arg(object):
    def __init__(self, typ, name):
        self.typ = str(typ)
        self.name = str(name)
    def type_fit(self, pol):
        # Modify the type to fit the owner.
        typ = self.typ
        if typ[-2] == "X":
            typ = typ[:-2] + pol[0] + typ[-1]
        if typ[-1] == "X":
            typ = typ[:-1] + pol[1]
        return typ
    def declare(self, pol):
        typ_fit = self.type_fit(pol)
        if typ_fit in POD_TYPES:
            return "{} {}".format(typ_fit, self.name)
        else:
            return "{} const& {}".format(typ_fit, self.name)
    def use(self):
        return self.name

class ConstructorCompose(object):
    def __init__(self, owner, child_pols):
        self.owner = owner
        self.child_pols = list(child_pols)
    def _args_declare(self):
        results = []
        for child_pol in self.child_pols:
            child_struct_name = struct_name(self.owner.base_name, child_pol)
            results.append("{} const& {}".format(
                child_struct_name,
                camel_to_snake(child_struct_name)))
        return ",\n".join(results)
    def declare(self):
        return "{}(\n{});\n".format(
            self.owner.name(),
            indent(self._args_declare()))
    def define(self):
        results = []
        for field in self.owner.fields:
            # Find which child also fits the mask.
            for child_pol in self.child_pols:
                child_struct_name = struct_name(self.owner.base_name, child_pol)
                if mask_match(child_pol, field.pol):
                    results.append("{0}({1}.{0})".format(
                        field.name,
                        camel_to_snake(child_struct_name)))
        return "inline {0}::{0}(\n{1}) :\n{2} {{ }}\n".format(
            self.owner.name(),
            indent(self._args_declare()),
            indent(",\n".join(results)))

class ConstructorDecompose(object):
    def __init__(self, owner, parent_pol):
        self.owner = owner
        self.parent_pol = list(parent_pol)
    def _args_declare(self):
        results = []
        parent_struct_name = struct_name(self.owner.base_name, self.parent_pol)
        return "{} const& {}".format(
            parent_struct_name,
            camel_to_snake(parent_struct_name))
    def declare(self):
        return "{}({});\n".format(self.owner.name(), self._args_declare())
    def define(self):
        results = []
        for field in self.owner.fields:
            parent_struct_name = struct_name(self.owner.base_name, self.parent_pol)
            results.append("{0}({1}.{0})".format(
                field.name,
                camel_to_snake(parent_struct_name)))
        return "inline {0}::{0}({1}) :\n{2} {{ }}\n".format(
            self.owner.name(),
            self._args_declare(),
            indent(",\n".join(results)))

class ConstructorDelegate(object):
    def __init__(self, owner, child_pols, args):
        self.owner = owner
        self.child_pols = list(child_pols)
        self.args = list(args)
    def _args_declare(self):
        results = []
        for arg in self.args:
            results.append(arg.declare(self.owner.pol))
        return ", ".join(results)
    def _args_use(self):
        results = []
        for arg in self.args:
            results.append(arg.use())
        return ", ".join(results)
    def declare(self):
        return "explicit {}({});\n".format(
            self.owner.name(),
            self._args_declare())
    def define(self):
        results = []
        pols = []
        for child_pol in self.child_pols:
            child_struct_name = struct_name(self.owner.base_name, child_pol)
            results.append("{}({})".format(
                child_struct_name,
                self._args_use()))
        return "inline {0}::{0}({1}) : {0}(\n{2}) {{ }}\n".format(
            self.owner.name(),
            self._args_declare(),
            indent(",\n".join(results)))

class ConstructorBase(object):
    def __init__(self, owner, args):
        self.owner = owner
        self.args = list(args)
    def _args_declare(self):
        results = []
        for arg in self.args:
            results.append(arg.declare(self.owner.pol))
        return ", ".join(results)
    def _args_use(self):
        results = []
        for arg in self.args:
            results.append(arg.use())
        return ", ".join(results)
    def declare(self):
        return "explicit {}({});\n".format(
            self.owner.name(),
            self._args_declare())
    def define(self):
        return ""

class ConstructorFields(object):
    def __init__(self, owner):
        self.owner = owner;
    def _args_declare(self):
        results = []
        for field in self.owner.fields:
            results.append(field.declare_as_arg())
        return ",\n".join(results)
    def declare(self):
        return "{}(\n{});\n".format(
            self.owner.name(),
           indent(self._args_declare()))
    def define(self):
        results = []
        for field in self.owner.fields:
            results.append("{0}({0})".format(field.use()))
        return "inline {0}::{0}(\n{1}) :\n{2} {{ }}\n".format(
            self.owner.name(),
            indent(self._args_declare()),
            indent(",\n".join(results)))

class ConstructorDefault(object):
    def __init__(self, owner):
        self.owner = owner;
    def declare(self):
        return "{}();\n".format(self.owner.name())
    def define(self):
        results = []
        for field in self.owner.fields:
            results.append("{0}()".format(field.use()))
        return "inline {0}::{0}() :\n{1} {{ }}\n".format(
            self.owner.name(),
            indent(",\n".join(results)))

class Struct(object):
    def __init__(self, base_name, pol, base_fields):
        self.base_name = str(base_name)
        self.pol = list(pol)
        self.fields = []
        for field_pol, fields in base_fields.items():
            if mask_match(self.pol, field_pol):
                for field in fields:
                    self.fields.append(Field(field_pol, field[0], field[1]))
        self.constructor_compose = []
        self.constructor_decompose = []
        self.constructor_delegate = []
        self.constructor_base = []
        self.constructor_fields = None
        self.constructor_default = None
    def name(self):
        return struct_name(self.base_name, self.pol)
    def declare(self):
        return "struct {};\n".format(self.name())
    def define(self):
        result = ""
        result += "struct {} {{\n".format(self.name())
        for field in self.fields:
            result += indent("{};\n".format(field.declare_as_field()))
        for cons in self.constructor_compose:
            result += indent(cons.declare())
        for cons in self.constructor_decompose:
            result += indent(cons.declare())
        for cons in self.constructor_delegate:
            result += indent(cons.declare())
        for cons in self.constructor_base:
            result += indent(cons.declare())
        if self.constructor_fields is not None:
            result += indent(self.constructor_fields.declare())
        if self.constructor_default is not None:
            result += indent(self.constructor_default.declare())
        result += "};\n"
        return result
    def define_methods(self):
        result = ""
        for cons in self.constructor_compose:
            result += cons.define()
        for cons in self.constructor_decompose:
            result += cons.define()
        for cons in self.constructor_delegate:
            result += cons.define()
        for cons in self.constructor_base:
            result += cons.define()
        if self.constructor_fields is not None:
            result += self.constructor_fields.define()
        if self.constructor_default is not None:
            result += self.constructor_default.define()
        return result
    def add_constructor_compose(self, child_pols):
        self.constructor_compose.append(
            ConstructorCompose(self, child_pols))
    def add_constructor_decompose(self, parent_pol):
        self.constructor_decompose.append(
            ConstructorDecompose(self, parent_pol))
    def add_constructor_delegate(self, child_pols, args):
        self.constructor_delegate.append(
            ConstructorDelegate(self, child_pols, args))
    def add_constructor_base(self, args):
        self.constructor_base.append(
            ConstructorBase(self, args))
    def add_constructor_fields(self):
        self.constructor_fields = ConstructorFields(self)
    def add_constructor_default(self):
        self.constructor_default = ConstructorDefault(self)

# Generates a set of structs for all combinations of beam and target
# polarizations, as well as conversions between them.
def generate_structs_pol(
        base_name,
        beam_pols,
        target_pols,
        base_fields,
        constructors=None,
        constructor_fields=None,
        constructor_default=None,
        generate_target_p=None):
    constructor_args = []
    if constructors is not None:
        for args in constructors:
            constructor_args.append([Arg(arg[0], arg[1]) for arg in args])
    structs = []
    # Structs with defined beam and target polarizations.
    for beam_pol in beam_pols:
        for target_pol in target_pols:
            pol = [beam_pol, target_pol]
            struct = Struct(base_name, pol, base_fields)
            if generate_target_p and target_pol != "U":
                struct.add_constructor_decompose([beam_pol, "P"])
            struct.add_constructor_decompose([beam_pol, "X"])
            struct.add_constructor_decompose(["X", target_pol])
            if generate_target_p and target_pol != "U":
                struct.add_constructor_decompose(["X", "P"])
            struct.add_constructor_decompose(["X", "X"])
            for args in constructor_args:
                struct.add_constructor_base(args)
            structs.append(struct)
    # Unknown target polarization structs.
    if generate_target_p:
        for beam_pol in beam_pols:
            pol = [beam_pol, "P"]
            struct = Struct(base_name, pol, base_fields)
            struct.add_constructor_compose(
                mask_find_all(pol, [beam_pol], target_pols))
            struct.add_constructor_decompose([beam_pol, "X"])
            struct.add_constructor_decompose(["X", "P"])
            struct.add_constructor_decompose(["X", "X"])
            for args in constructor_args:
                struct.add_constructor_delegate(
                    mask_find_all(pol, [beam_pol], target_pols),
                    args)
            structs.append(struct)
    for beam_pol in beam_pols:
        pol = [beam_pol, "X"]
        struct = Struct(base_name, pol, base_fields)
        struct.add_constructor_compose(
            mask_find_all(pol, [beam_pol], target_pols))
        if generate_target_p:
            struct.add_constructor_compose([[beam_pol, "U"], [beam_pol, "P"]])
        struct.add_constructor_decompose(["X", "X"])
        for args in constructor_args:
            struct.add_constructor_delegate(
                mask_find_all(pol, [beam_pol], target_pols),
                args)
        structs.append(struct)
    # Unknown beam polarization structs.
    for target_pol in target_pols:
        pol = ["X", target_pol]
        struct = Struct(base_name, pol, base_fields)
        struct.add_constructor_compose(
            mask_find_all(pol, beam_pols, [target_pol]))
        if generate_target_p and target_pol != "U":
            struct.add_constructor_decompose(["X", "P"])
        struct.add_constructor_decompose(["X", "X"])
        for args in constructor_args:
            struct.add_constructor_delegate(
                mask_find_all(pol, beam_pols, [target_pol]),
                args)
        structs.append(struct)
    # Generally polarized structs.
    if generate_target_p:
        struct = Struct(base_name, ["X", "P"], base_fields)
        struct.add_constructor_compose(
            mask_find_all(["X", "P"], beam_pols, target_pols))
        struct.add_constructor_compose(
            mask_find_all(["X", "P"], beam_pols, ["P"]))
        struct.add_constructor_compose(
            mask_find_all(["X", "P"], ["X"], target_pols))
        struct.add_constructor_decompose(["X", "X"])
        for args in constructor_args:
            struct.add_constructor_delegate(
                mask_find_all(["X", "P"], beam_pols, target_pols),
                args)
        structs.append(struct)
    struct = Struct(base_name, ["X", "X"], base_fields)
    struct.add_constructor_compose(
        mask_find_all(["X", "X"], beam_pols, target_pols))
    if generate_target_p:
        struct.add_constructor_compose(
            mask_find_all(["X", "X"], beam_pols, ["U", "P"]))
    struct.add_constructor_compose(
        mask_find_all(["X", "X"], beam_pols, ["X"]))
    if generate_target_p:
        struct.add_constructor_compose(
            mask_find_all(["X", "X"], ["X"], ["U", "P"]))
    struct.add_constructor_compose(
        mask_find_all(["X", "X"], ["X"], target_pols))
    for args in constructor_args:
        struct.add_constructor_delegate(
            mask_find_all(["X", "X"], beam_pols, target_pols),
            args)
    structs.append(struct)

    # Add optional field constructors.
    if constructor_fields:
        for struct in structs:
            struct.add_constructor_fields()
    if constructor_default:
        for struct in structs:
            struct.add_constructor_default()

    # Output the structs.
    for struct in structs:
        cog.out(struct.declare())
    cog.out("\n")
    for struct in structs:
        cog.out(struct.define())
        cog.out("\n")
    for struct in structs:
        cog.out(struct.define_methods())
        cog.out("\n")

    # Type aliases.
    cog.out("using {0} = {0}XX;\n".format(base_name))

