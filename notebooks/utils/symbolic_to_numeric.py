import re
import sympy as sp
from itertools import product

import os
this_dir = os.path.dirname(__file__)  # carpeta donde está este .py
template_path_metric = os.path.join(this_dir, "metric_template.txt")
template_path_coordinates = os.path.join(this_dir, "coordinate_template.txt")

template_metric = open(template_path_metric).read()
template_coordinates = open(template_path_coordinates).read()


def make_valid_varname(name):
    # reemplazar caracteres no válidos por "_"
    # {'**':'pow', '(':l, ')':l}
    replaces = {
        r"\*\*": "_pow",
        r"\(": "_I",
        r"\)": "I_",
        r"/": "_over_",
        r"\+": "_plus_"
    }

    # replace '\dot{x}^i' by 'dxi_dt'
    name = re.sub(r'\\dot\{x\}\^(-?\d+)', r'dx\1_dt', name)

    for pattern, replacement in replaces.items():
        name = re.sub(pattern, replacement, name)

    # _powN1_powN2 -> _pow#(N1*N2)
    name = re.sub(r'_pow(-?\d+)_pow(-?\d+)', lambda m: f"_pow{int(m.group(1)) * int(m.group(2))}", name)
    name = re.sub(r'[^0-9a-zA-Z_]', '_', name)
    # no empezar con dígito
    if name[0].isdigit():
        name = "_" + name
    return name

def generate_new_symbols(original_symbols):
    new_symbols = {}
    for symbol in original_symbols:
    # reemplazamos 'x^i' por 'xi'
        str_symbol = str(symbol)
        if 'x^' in str_symbol:
            str_symbol = str_symbol.replace('x^', 'x')
        else:
            str_symbol = make_valid_varname(str_symbol)

        new_symbols[symbol] = sp.symbols(str_symbol)

    return new_symbols


def rename_symbols_of_expr(expr):
    original_symbols = expr.free_symbols
    coordinate_substitutions = generate_new_symbols(original_symbols)
    new_expr = expr.subs(coordinate_substitutions)

    return new_expr, coordinate_substitutions

def get_optimized_christoffel_symbols(metric):
    """
    Optimiza los símbolos de Christoffel mediante sustituciones simbólicas
    y eliminación de componentes simétricas redundantes.
    
    Args:
        metric: Tensor de métrica simbólica

    Returns:
        tuple: (símbolos optimizados, lista de sustituciones aplicadas)
    """
    christoffel_symbols = metric.christoffel_symbols().tensor()

    substitution_history = []

    # Primera optimización: sustituir símbolos originales 
    # original_symbols = christoffel_symbols.free_symbols
    # coordinate_substitutions = generate_new_symbols(original_symbols)
    
    # optimized_symbols = christoffel_symbols.subs(coordinate_substitutions)
    optimized_symbols, coordinate_substitutions = rename_symbols_of_expr(christoffel_symbols)

    # Guardar sustituciones excluyendo coordenadas x^i
    non_coordinate_subs = {
        symbol: replacement 
        for symbol, replacement in coordinate_substitutions.items() 
        if 'x^' not in str(symbol)
    }
    substitution_history.append(non_coordinate_subs)

    # Convertir a array mutable para manipulación eficiente
    symbol_array = sp.MutableDenseNDimArray(
        optimized_symbols.tolist(), 
        (4, 4, 4)
    )

    # Eliminar componentes simétricas redundantes (mantener solo nu <= sigma)
    for mu, nu, sigma in product(range(4), repeat=3):
        if nu > sigma:
            symbol_array[mu, nu, sigma] = 0

    # Restaurar tipo original del tensor
    optimized_symbols = optimized_symbols.__class__(
        symbol_array.tolist(), 
        (4, 4, 4)
    )

    # Optimizaciones por tipo de expresión matemática
    expression_types = [sp.Function, sp.Pow, sp.Pow, sp.Mul]

    for expr_type in expression_types:
        # Encontrar expresiones repetidas del tipo actual
        repeated_expressions = _find_repeated_expressions(optimized_symbols, expr_type)
        
        if repeated_expressions:
            # Generar nuevas variables simbólicas para las expresiones repetidas
            expression_substitutions = generate_new_symbols(repeated_expressions)
            optimized_symbols = optimized_symbols.subs(expression_substitutions)
            substitution_history.append(expression_substitutions)

    return optimized_symbols, substitution_history

def _find_repeated_expressions(tensor, expression_type):
    """
    Encuentra expresiones del tipo especificado que aparecen múltiples veces
    y no contienen sumas en sus argumentos.
    
    Args:
        tensor: Tensor a analizar
        expression_type: Tipo de expresión a buscar (sp.Function, sp.Pow, etc.)
        
    Returns:
        set: Conjunto de expresiones repetidas
    """
    # Extraer todas las expresiones del tipo especificado
    all_expressions = {
        expr for expr in tensor.atoms(expression_type)
        # if not any(isinstance(arg, sp.Add) for arg in expr.args)
    }
    
    # Filtrar solo las que aparecen más de una vez
    repeated_expressions = {
        expr for expr in all_expressions 
        if tensor.count(expr) > 1
    }
    
    return repeated_expressions


def fill_template(template: str, replacements: dict) -> str:
    """
    Reemplaza placeholders en un template con los valores provistos.
    Los placeholders deben estar en la forma $key$ dentro del template.
    La indentación del placeholder se aplicará a todas las líneas del contenido de reemplazo.

    :param template: Texto base con placeholders.
    :param replacements: Diccionario { "key": "contenido a insertar" }.
    :return: Texto expandido.
    """
    def replace(match):
        # Capturar la indentación antes del placeholder
        full_match = match.group(0)
        key = match.group(2)
        leading_spaces = match.group(1)[1:]
        
        if key not in replacements:
            return full_match  # si no existe, dejar igual
        
        content = replacements[key]
        
        # Si el contenido está vacío, retornar vacío
        if not content:
            return ""
        
        # Dividir el contenido en líneas
        lines = content.split('\n')
        
        # Si solo hay una línea, no necesitamos procesar indentación adicional
        if len(lines) == 1:
            return content
        
        # Aplicar la indentación del placeholder a todas las líneas (incluida la primera)
        indented_lines = []
        
        for line in lines:
            if line.strip():  # Si la línea no está vacía
                indented_lines.append(leading_spaces + line)
            else:  # Línea vacía, mantener vacía
                indented_lines.append('')
        
        return '\n'.join(indented_lines)
    
    # Capturar espacios/tabs antes del placeholder
    return re.sub(r'^(\s*)\$(\w+)\$$', replace, template, flags=re.MULTILINE)


def export_christoffel_symbols_to_code(metric):
    """
    Exporta los símbolos de Christoffel optimizados a código Python.
    
    Args:
        ch_subs: Tensor de símbolos de Christoffel optimizados
        all_subs: Lista de diccionarios con las sustituciones aplicadas
    
    Returns:
        str: Código Python con las definiciones de variables y asignaciones
    """
    ch_subs, all_subs = get_optimized_christoffel_symbols(metric)

    # Ahora puedes construir tus variables sin preocuparte por la indentación
    auxiliary_variables = ""

    mega_dict = all_subs[0].copy()

    for key, value in all_subs[0].items():
        auxiliary_variables += f"{value} = self.{key}\n"
    auxiliary_variables += "\n"

    for i, sub_ in enumerate(all_subs[1:]):
        sub_ = dict(sorted(sub_.items(), key=lambda item: str(item[1])))
        for key, value in sub_.items():
            if value not in mega_dict.values():
                auxiliary_variables += f"{value} = {key}\n"
        i += 1
        mega_dict.update(sub_)
        auxiliary_variables += "\n"

    auxiliary_variables += "\n"

    gamma_txt = ""

    for mu in range(0, 4):
        gamma_txt += "\n"
        gamma_txt += f"# mu = {mu}\n"
        for nu, rho in product(range(4), repeat=2):
            if nu <= rho:
                if ch_subs[mu, nu, rho] != 0:
                    gamma_txt += f"Gamma[{mu}, {nu}, {rho}] = {ch_subs[mu, nu, rho]}\n"
        gamma_txt += "\n"
        for nu, rho in product(range(4), repeat=2):
            if nu > rho:
                if ch_subs[mu, rho, nu] != 0:
                    gamma_txt += f"Gamma[{mu}, {nu}, {rho}] = Gamma[{mu}, {rho}, {nu}]\n"

    # Usar el template
    txt = fill_template(template_metric, {
        'auxiliary_variables': auxiliary_variables,
        'gamma_txt': gamma_txt
    })

    return txt