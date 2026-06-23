from models import Reaction

def format_reaction(reaction: Reaction) -> str:
    reactants_list = []
    products_list = []
    
    for coeff, value in reaction.reactants: 
        term = f"{coeff} {value}" if coeff != 1 else str(value)
        reactants_list.append(term)
                
            
    for coeff, value in reaction.products: 
        term = f"{coeff} {value}" if coeff != 1 else str(value)
        products_list.append(term)

    return " + ".join(reactants_list) + " ---> " + " + ".join(products_list)

def format_reaction_list(reactions: list[Reaction]) -> str:
    lines = []
    for i, reaction in enumerate(reactions, 1):
        lines.append(f"{i}. {format_reaction(reaction)}")
    return "\n".join(lines)
