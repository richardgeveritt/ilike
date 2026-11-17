#ifndef ILIKE_LANGUAGE_H
#define ILIKE_LANGUAGE_H

// The parser of the ilike language (dev/ilike-language.md; dev/plan.md 3.14,
// I3a). It reads a file into a syntax tree, as a RecipeValue, without
// interpreting the keywords: the sections, their statements
// ([keyword [type arguments] [<constraints>]] [name] (= | ~) expression),
// expressions, and the code sections verbatim. Interpreting the tree (the
// recipe the engine reads) is the next step (I3b). Errors give
// file:line:column.

#include <string>
#include "recipe_value.h"

namespace ilike
{

// The syntax tree of one file:
// {version, sections: [section]}, a section being
//   {section: "model" | "constants" | "data" | "parameters" | "algorithm",
//    name (algorithm), statements: [statement], line, column}
//   {section: "cpp" | "r" | "py" | "jl", code, line, column}
//   {section: "import", path, alias, line, column}
// a statement
//   {keyword, type_arguments: [expression], constraints: [named expressions],
//    target: expression (a name, or a name with indices), operator ("=", "~"
//    "<-" for a deterministic node, or "" for a declaration), value:
//    expression, line, column}
//   {keyword: "for", index, range: expression, body: [statement], line, column}
// and an expression
//   {expression: "number", value}  {expression: "string", value}
//   {expression: "logical", value}  {expression: "name", name, namespace}
//   {expression: "call", function, namespace, arguments (named: keyword)}
//   {expression: "vector", elements}  {expression: "index", target, indices}
//   {expression: "binary", operator, left, right}
//   {expression: "unary", operator, operand}
//   {expression: "member", target, name}   (target$name, for other languages)
// each with line and column.
RecipeValue parse_ilike_language(const std::string &text, const std::string &file_name);

// Lowering (I3b): from the syntax trees of the files of a recipe (a list of
// trees) to {recipe (the model: data and factors), algorithms (a named list,
// each with its moves, settings of the old method blocks, and run-time
// settings), constants (defaults by name), user_code (the cpp sections),
// generated_code (small functions made by the lowering)}. User functions are
// placeholders {ilike_function: name, role: old block name}.
RecipeValue lower_ilike_language(const RecipeValue &trees);

// The signature "type name(arguments)" of the definition of the C++ function
// name in code, or "" if it is not defined there.
std::string find_cpp_function_signature(const std::string &code, const std::string &name);

// TRUE if the text is in the ilike language (its first statement, after
// comments, is the version line "ilike <version>").
bool is_ilike_language(const std::string &text);

}

#endif
