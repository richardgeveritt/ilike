// The parser of the ilike language (dev/ilike-language.md; dev/plan.md 3.14,
// I3a): a hand-written lexer and recursive-descent parser producing a syntax
// tree (RecipeValue). See ilike_language.h for the tree's form.

#include <RcppArmadillo.h>

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <deque>
#include <sstream>

#include "ilike_language.h"

namespace ilike
{

namespace
{

enum class TokenType { Identifier, Number, String, Symbol, Newline, End };

struct Token
{
  TokenType type = TokenType::End;
  std::string text;
  double number = 0.0;
  size_t line = 1;
  size_t column = 1;
};

// The lexer reads tokens on demand, so that a code section can be taken
// verbatim from the text (capture_code) once its opening brace is reached.
// Newlines are tokens only outside parentheses and brackets (so expressions
// may continue over several lines inside them).
class Lexer
{

public:

  Lexer(const std::string &text_in, const std::string &file_in)
    : text(text_in), file(file_in) {}

  const Token &peek(size_t k = 0)
  {
    while (this->buffer.size() <= k)
      this->buffer.push_back(this->read());
    return this->buffer[k];
  }

  Token next()
  {
    Token token = this->peek();
    this->buffer.pop_front();
    return token;
  }

  [[noreturn]] void error(size_t line, size_t column, const std::string &message) const
  {
    Rcpp::stop(this->file + ":" + std::to_string(line) + ":" + std::to_string(column) + ": " + message);
  }

  // After the opening brace of a code section (with nothing read ahead):
  // the text up to the matching closing brace, respecting the language's
  // strings and comments.
  std::string capture_code(const std::string &language, size_t start_line, size_t start_column)
  {
    if (!this->buffer.empty())
      this->error(start_line, start_column, "internal error: tokens read ahead of a code section.");
    size_t depth = 1;
    size_t start = this->position;
    bool c_like = (language == "cpp");
    while (this->position < this->text.size())
    {
      char c = this->current();
      char d = this->ahead(1);
      if (c_like && c == '/' && d == '/')
      {
        while (this->position < this->text.size() && this->current() != '\n')
          this->advance();
        continue;
      }
      if (c_like && c == '/' && d == '*')
      {
        this->advance(); this->advance();
        while (this->position < this->text.size() && !(this->current() == '*' && this->ahead(1) == '/'))
          this->advance();
        if (this->position >= this->text.size())
          this->error(start_line, start_column, "unterminated comment in the " + language + " section.");
        this->advance(); this->advance();
        continue;
      }
      if (!c_like && c == '#')
      {
        if (language == "jl" && d == '=')
        {
          this->advance(); this->advance();
          while (this->position < this->text.size() && !(this->current() == '=' && this->ahead(1) == '#'))
            this->advance();
          if (this->position >= this->text.size())
            this->error(start_line, start_column, "unterminated comment in the " + language + " section.");
          this->advance(); this->advance();
          continue;
        }
        while (this->position < this->text.size() && this->current() != '\n')
          this->advance();
        continue;
      }
      if (c_like && c == 'R' && d == '"' &&
          (this->position == 0 || !(std::isalnum((unsigned char)this->text[this->position - 1]) ||
                                    this->text[this->position - 1] == '_')))
      {
        // a C++ raw string R"delimiter( ... )delimiter"
        size_t open = this->text.find('(', this->position + 2);
        if (open != std::string::npos)
        {
          std::string delimiter = ")" + this->text.substr(this->position + 2, open - this->position - 2) + "\"";
          size_t close = this->text.find(delimiter, open + 1);
          if (close == std::string::npos)
            this->error(this->line, this->column, "unterminated raw string in the " + language + " section.");
          while (this->position < close + delimiter.size())
            this->advance();
          continue;
        }
      }
      if (c == '"' || c == '\'')
      {
        // (in Python, triple-quoted strings)
        bool triple = (language == "py") && (d == c) && (this->ahead(2) == c);
        size_t quote_line = this->line, quote_column = this->column;
        this->advance();
        if (triple) { this->advance(); this->advance(); }
        while (true)
        {
          if (this->position >= this->text.size())
            this->error(quote_line, quote_column, "unterminated string in the " + language + " section.");
          char e = this->current();
          if (e == '\\') { this->advance(); this->advance(); continue; }
          if (triple && e == c && this->ahead(1) == c && this->ahead(2) == c)
          {
            this->advance(); this->advance(); this->advance();
            break;
          }
          if (!triple && e == c) { this->advance(); break; }
          if (!triple && e == '\n' && c_like)
            this->error(quote_line, quote_column, "unterminated string in the " + language + " section.");
          this->advance();
        }
        continue;
      }
      if (c == '{')
        ++depth;
      else if (c == '}')
      {
        --depth;
        if (depth == 0)
        {
          std::string code = this->text.substr(start, this->position - start);
          this->advance();
          return code;
        }
      }
      this->advance();
    }
    this->error(start_line, start_column, "the " + language + " section has no closing brace.");
  }

  // the raw character after the next token's whitespace and comments
  // (newlines included), with nothing read ahead
  char next_raw_character()
  {
    this->skip_space(true);
    return (this->position < this->text.size()) ? this->current() : '\0';
  }

  void consume_raw_character()
  {
    this->advance();
  }

  size_t current_line() const { return this->line; }
  size_t current_column() const { return this->column; }

private:

  std::string text;
  std::string file;
  size_t position = 0;
  size_t line = 1;
  size_t column = 1;
  int depth = 0;               // parentheses and brackets
  std::deque<Token> buffer;

  char current() const { return this->text[this->position]; }
  char ahead(size_t k) const
  {
    return (this->position + k < this->text.size()) ? this->text[this->position + k] : '\0';
  }

  void advance()
  {
    if (this->position >= this->text.size())
      return;
    if (this->text[this->position] == '\n')
    {
      ++this->line;
      this->column = 1;
    }
    else
      ++this->column;
    ++this->position;
  }

  // spaces and comments; newlines too if include_newlines (or inside
  // parentheses and brackets)
  void skip_space(bool include_newlines)
  {
    while (this->position < this->text.size())
    {
      char c = this->current();
      if (c == '\n' && !include_newlines && this->depth == 0)
        return;
      if (std::isspace((unsigned char)c))
      {
        this->advance();
        continue;
      }
      if (c == '/' && this->ahead(1) == '/')
      {
        while (this->position < this->text.size() && this->current() != '\n')
          this->advance();
        continue;
      }
      if (c == '/' && this->ahead(1) == '*')
      {
        size_t comment_line = this->line, comment_column = this->column;
        this->advance(); this->advance();
        while (this->position < this->text.size() && !(this->current() == '*' && this->ahead(1) == '/'))
          this->advance();
        if (this->position >= this->text.size())
          this->error(comment_line, comment_column, "unterminated comment.");
        this->advance(); this->advance();
        continue;
      }
      return;
    }
  }

  Token read()
  {
    this->skip_space(false);
    Token token;
    token.line = this->line;
    token.column = this->column;
    if (this->position >= this->text.size())
    {
      token.type = TokenType::End;
      return token;
    }
    char c = this->current();
    if (c == '\n')
    {
      token.type = TokenType::Newline;
      while (this->position < this->text.size())
      {
        this->skip_space(false);
        if (this->position < this->text.size() && this->current() == '\n')
          this->advance();
        else
          break;
      }
      return token;
    }
    if (std::isalpha((unsigned char)c) || c == '_' || (unsigned char)c >= 0x80)
    {
      // (non-ASCII bytes: UTF-8 letters such as θ, as R allows in names)
      size_t start = this->position;
      while (this->position < this->text.size() &&
             (std::isalnum((unsigned char)this->current()) || this->current() == '_' ||
              (unsigned char)this->current() >= 0x80))
        this->advance();
      token.type = TokenType::Identifier;
      token.text = this->text.substr(start, this->position - start);
      return token;
    }
    if (std::isdigit((unsigned char)c) || (c == '.' && std::isdigit((unsigned char)this->ahead(1))))
    {
      const char *begin = this->text.c_str() + this->position;
      char *end = nullptr;
      token.number = std::strtod(begin, &end);
      size_t length = end - begin;
      token.type = TokenType::Number;
      token.text = this->text.substr(this->position, length);
      for (size_t k = 0; k < length; ++k)
        this->advance();
      // R's integer suffix: 3L (marked, for calls to other languages)
      if (this->position < this->text.size() && this->current() == 'L' &&
          token.text.find_first_of(".eE") == std::string::npos &&
          !(this->position + 1 < this->text.size() &&
            (std::isalnum((unsigned char)this->text[this->position + 1]) || this->text[this->position + 1] == '_')))
      {
        token.text += "L";
        this->advance();
      }
      if (this->position < this->text.size() &&
          (std::isalpha((unsigned char)this->current()) || this->current() == '_'))
        this->error(token.line, token.column, "malformed number \"" + token.text + this->current() + "\".");
      return token;
    }
    if (c == '"')
    {
      this->advance();
      std::string value;
      while (true)
      {
        if (this->position >= this->text.size() || this->current() == '\n')
          this->error(token.line, token.column, "unterminated string.");
        char e = this->current();
        if (e == '\\')
        {
          this->advance();
          char f = this->current();
          value.push_back(f == 'n' ? '\n' : (f == 't' ? '\t' : f));
          this->advance();
          continue;
        }
        if (e == '"')
        {
          this->advance();
          break;
        }
        value.push_back(e);
        this->advance();
      }
      token.type = TokenType::String;
      token.text = value;
      return token;
    }
    // symbols
    token.type = TokenType::Symbol;
    if (c == ':' && this->ahead(1) == ':')
    {
      token.text = "::";
      this->advance(); this->advance();
      return token;
    }
    static const std::string singles = "()[]{},;=~+-*/^<>:$";
    if (singles.find(c) == std::string::npos)
      this->error(token.line, token.column, std::string("unexpected character '") + c + "'.");
    token.text = std::string(1, c);
    if (c == '(' || c == '[')
      ++this->depth;
    else if ((c == ')' || c == ']') && this->depth > 0)
      --this->depth;
    this->advance();
    return token;
  }

};

// --- the parser ----------------------------------------------------------------

class Parser
{

public:

  Parser(const std::string &text, const std::string &file_in)
    : lexer(text, file_in), file(file_in) {}

  RecipeValue program()
  {
    this->skip_newlines();
    Token version_word = this->lexer.next();
    if (!(version_word.type == TokenType::Identifier && version_word.text == "ilike"))
      this->lexer.error(version_word.line, version_word.column,
                        "an ilike-language file begins with the version line \"ilike 1.0\".");
    Token version = this->lexer.next();
    if (version.type != TokenType::Number)
      this->lexer.error(version.line, version.column, "expected the language version after \"ilike\" (e.g. ilike 1.0).");
    if (std::floor(version.number) != 1.0)
      this->lexer.error(version.line, version.column, "this version of the package reads version 1 of the ilike language, not " + version.text + ".");
    this->end_of_statement();

    RecipeValue sections = RecipeValue::list();
    while (true)
    {
      this->skip_newlines();
      if (this->lexer.peek().type == TokenType::End)
        break;
      sections.push_back(this->section());
    }

    RecipeValue output = RecipeValue::list();
    output.push_back(RecipeValue(version.text), "version");
    output.push_back(RecipeValue(this->file), "file");
    output.push_back(sections, "sections");
    return output;
  }

private:

  Lexer lexer;
  std::string file;

  static RecipeValue position(const Token &token)
  {
    RecipeValue output = RecipeValue::list();
    output.push_back(RecipeValue(int(token.line)), "line");
    output.push_back(RecipeValue(int(token.column)), "column");
    return output;
  }

  static void add_position(RecipeValue &node, const Token &token)
  {
    node.push_back(RecipeValue(int(token.line)), "line");
    node.push_back(RecipeValue(int(token.column)), "column");
  }

  static std::string describe(const Token &token)
  {
    switch (token.type)
    {
      case TokenType::Identifier: return "\"" + token.text + "\"";
      case TokenType::Number: return token.text;
      case TokenType::String: return "the string \"" + token.text + "\"";
      case TokenType::Symbol: return "\"" + token.text + "\"";
      case TokenType::Newline: return "the end of the line";
      case TokenType::End: return "the end of the file";
    }
    return "";
  }

  bool is_symbol(const Token &token, const std::string &symbol) const
  {
    return token.type == TokenType::Symbol && token.text == symbol;
  }

  Token expect_symbol(const std::string &symbol, const std::string &context)
  {
    Token token = this->lexer.next();
    if (!this->is_symbol(token, symbol))
      this->lexer.error(token.line, token.column, "expected \"" + symbol + "\" " + context + ", found " + describe(token) + ".");
    return token;
  }

  Token expect_identifier(const std::string &context)
  {
    Token token = this->lexer.next();
    if (token.type != TokenType::Identifier)
      this->lexer.error(token.line, token.column, "expected " + context + ", found " + describe(token) + ".");
    return token;
  }

  void skip_newlines()
  {
    while (this->lexer.peek().type == TokenType::Newline ||
           this->is_symbol(this->lexer.peek(), ";"))
      this->lexer.next();
  }

  // a statement ends at a newline, ";", "}" (not consumed) or the end
  void end_of_statement()
  {
    const Token &token = this->lexer.peek();
    if (token.type == TokenType::Newline || this->is_symbol(token, ";"))
    {
      this->lexer.next();
      return;
    }
    if (token.type == TokenType::End || this->is_symbol(token, "}"))
      return;
    this->lexer.error(token.line, token.column, "expected the end of the statement, found " + describe(token) + ".");
  }

  RecipeValue section()
  {
    Token keyword = this->expect_identifier("a section (constants, data, parameters, model, algorithm, cpp, r, py, jl or import)");
    const std::string &kind = keyword.text;
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue(kind), "section");

    if (kind == "import")
    {
      Token path = this->lexer.next();
      if (path.type != TokenType::String)
        this->lexer.error(path.line, path.column, "expected the file to import, in quotes.");
      Token as = this->expect_identifier("\"as\" after the file to import");
      if (as.text != "as")
        this->lexer.error(as.line, as.column, "expected \"as\" after the file to import.");
      Token alias = this->expect_identifier("the name to import the file as");
      node.push_back(RecipeValue(path.text), "path");
      node.push_back(RecipeValue(alias.text), "alias");
      add_position(node, keyword);
      this->end_of_statement();
      return node;
    }

    if (kind == "cpp" || kind == "r" || kind == "py" || kind == "jl")
    {
      // nothing may be read ahead before the code is taken verbatim
      if (this->lexer.next_raw_character() != '{')
        this->lexer.error(this->lexer.current_line(), this->lexer.current_column(),
                          "expected \"{\" to begin the " + kind + " section.");
      size_t line = this->lexer.current_line(), column = this->lexer.current_column();
      this->lexer.consume_raw_character();
      std::string code = this->lexer.capture_code(kind, line, column);
      node.push_back(RecipeValue(code), "code");
      add_position(node, keyword);
      this->end_of_statement();
      return node;
    }

    if (kind == "algorithm")
    {
      std::string name;
      if (this->lexer.peek().type == TokenType::Identifier)
        name = this->lexer.next().text;
      node.push_back(RecipeValue(name), "name");
    }
    else if (!(kind == "constants" || kind == "data" || kind == "parameters" || kind == "model"))
    {
      this->lexer.error(keyword.line, keyword.column,
                        "unknown section \"" + kind + "\" (sections: constants, data, parameters, model, algorithm, cpp, r, py, jl, import).");
    }

    this->expect_symbol("{", "to begin the " + kind + " section");
    RecipeValue statements = RecipeValue::list();
    while (true)
    {
      this->skip_newlines();
      const Token &token = this->lexer.peek();
      if (this->is_symbol(token, "}"))
      {
        this->lexer.next();
        break;
      }
      if (token.type == TokenType::End)
        this->lexer.error(keyword.line, keyword.column, "the " + kind + " section has no closing brace.");
      statements.push_back(this->statement());
      this->end_of_statement();
    }
    node.push_back(statements, "statements");
    add_position(node, keyword);
    this->end_of_statement();
    return node;
  }

  // [keyword [ "[" arguments "]" ] [ "<" constraints ">" ]] target [(= | ~) value]
  RecipeValue statement()
  {
    Token first = this->lexer.peek();
    if (first.type != TokenType::Identifier)
      this->lexer.error(first.line, first.column, "expected a statement, found " + describe(first) + ".");
    if (first.text == "for")
      return this->for_statement();

    RecipeValue node = RecipeValue::list();
    std::string keyword;
    RecipeValue type_arguments = RecipeValue::list();
    RecipeValue constraints = RecipeValue::list();
    RecipeValue target;

    // the first word, then decide whether it is a keyword (a type or a kind
    // of statement) or the target
    Token word = this->lexer.next();
    RecipeValue word_indices;
    bool word_indexed = false;
    const Token *after = &this->lexer.peek();
    bool declaration_like = false;

    bool deterministic_next = this->is_symbol(*after, "<") && this->is_symbol(this->lexer.peek(1), "-") &&
                              this->lexer.peek(1).line == after->line && this->lexer.peek(1).column == after->column + 1;
    if (this->is_symbol(*after, "[") || (this->is_symbol(*after, "<") && !deterministic_next))
    {
      // vector[2] x, real<lower=0> x, or an indexed target y[i] ~ ...
      while (true)
      {
        after = &this->lexer.peek();
        if (this->is_symbol(*after, "["))
        {
          this->lexer.next();
          word_indices = this->arguments("]");
          word_indexed = true;
          continue;
        }
        if (this->is_symbol(*after, "<") &&
            !(this->is_symbol(this->lexer.peek(1), "-") && this->lexer.peek(1).line == after->line &&
              this->lexer.peek(1).column == after->column + 1))  // (not "<-")
        {
          this->lexer.next();
          constraints = this->constraint_list();
          declaration_like = true;
          continue;
        }
        break;
      }
      after = &this->lexer.peek();
      if (after->type == TokenType::Identifier)
        declaration_like = true;
    }
    else if (after->type == TokenType::Identifier)
    {
      declaration_like = true;
    }

    if (declaration_like)
    {
      keyword = word.text;
      if (word_indexed)
        type_arguments = word_indices;
      Token name = this->expect_identifier("a name after \"" + keyword + "\"");
      target = this->name_node(name, "");
      // an indexed name in a declaration (Stan's "vector[n] x" is above;
      // "real x[n]" is accepted too)
      if (this->is_symbol(this->lexer.peek(), "["))
      {
        this->lexer.next();
        RecipeValue indices = this->arguments("]");
        RecipeValue indexed = RecipeValue::list();
        indexed.push_back(RecipeValue("index"), "expression");
        indexed.push_back(target, "target");
        indexed.push_back(indices, "indices");
        add_position(indexed, name);
        target = indexed;
      }
    }
    else
    {
      target = this->name_node(word, "");
      if (word_indexed)
      {
        RecipeValue indexed = RecipeValue::list();
        indexed.push_back(RecipeValue("index"), "expression");
        indexed.push_back(target, "target");
        indexed.push_back(word_indices, "indices");
        add_position(indexed, word);
        target = indexed;
      }
    }

    std::string op;
    RecipeValue value;
    const Token &next = this->lexer.peek();
    if (this->is_symbol(next, "=") || this->is_symbol(next, "~"))
    {
      op = this->lexer.next().text;
      this->skip_newlines_inside();
      value = this->expression();
    }
    else if (this->is_symbol(next, "<") && this->is_symbol(this->lexer.peek(1), "-") &&
             this->lexer.peek(1).line == next.line && this->lexer.peek(1).column == next.column + 1)
    {
      // a deterministic node: mu[i] <- expression
      this->lexer.next();
      this->lexer.next();
      op = "<-";
      this->skip_newlines_inside();
      value = this->expression();
    }
    else if (!declaration_like)
    {
      this->lexer.error(next.line, next.column, "expected \"=\" or \"~\" after \"" + word.text + "\", found " + describe(next) + ".");
    }

    node.push_back(RecipeValue(keyword), "keyword");
    node.push_back(type_arguments, "type_arguments");
    node.push_back(constraints, "constraints");
    node.push_back(target, "target");
    node.push_back(RecipeValue(op), "operator");
    node.push_back(value, "value");
    add_position(node, first);
    return node;
  }

  // newlines after "=", "~" or an operator continue the statement
  void skip_newlines_inside()
  {
    while (this->lexer.peek().type == TokenType::Newline)
      this->lexer.next();
  }

  // < name = expression {, name = expression} >
  RecipeValue constraint_list()
  {
    RecipeValue output = RecipeValue::list();
    while (true)
    {
      Token name = this->expect_identifier("a constraint (e.g. lower = 0)");
      this->expect_symbol("=", "after the constraint \"" + name.text + "\"");
      output.push_back(this->additive(), name.text);
      Token token = this->lexer.next();
      if (this->is_symbol(token, ">"))
        break;
      if (!this->is_symbol(token, ","))
        this->lexer.error(token.line, token.column, "expected \",\" or \">\" in the constraints, found " + describe(token) + ".");
    }
    return output;
  }

  // for (i in a:b) statement, or for (i in a:b) { statements }
  RecipeValue for_statement()
  {
    Token keyword = this->lexer.next();
    this->expect_symbol("(", "after \"for\"");
    Token index = this->expect_identifier("the loop's index after \"for (\"");
    Token in = this->expect_identifier("\"in\" after the loop's index");
    if (in.text != "in")
      this->lexer.error(in.line, in.column, "expected \"in\" after the loop's index, found " + describe(in) + ".");
    RecipeValue range = this->expression();
    this->expect_symbol(")", "to close the loop's range");
    this->skip_newlines_inside();
    RecipeValue body = RecipeValue::list();
    if (this->is_symbol(this->lexer.peek(), "{"))
    {
      Token open = this->lexer.next();
      while (true)
      {
        this->skip_newlines();
        const Token &token = this->lexer.peek();
        if (this->is_symbol(token, "}"))
        {
          this->lexer.next();
          break;
        }
        if (token.type == TokenType::End)
          this->lexer.error(open.line, open.column, "the loop's block has no closing brace.");
        body.push_back(this->statement());
        this->end_of_statement();
      }
    }
    else
      body.push_back(this->statement());
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue("for"), "keyword");
    node.push_back(RecipeValue(index.text), "index");
    node.push_back(range, "range");
    node.push_back(body, "body");
    add_position(node, keyword);
    return node;
  }

  RecipeValue name_node(const Token &token, const std::string &name_space)
  {
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue("name"), "expression");
    node.push_back(RecipeValue(token.text), "name");
    node.push_back(RecipeValue(name_space), "namespace");
    add_position(node, token);
    return node;
  }

  // expressions, by precedence: + -, * /, : (ranges), unary -, ^, postfix
  RecipeValue expression() { return this->additive(); }

  RecipeValue binary(const std::string &op, const RecipeValue &left, const RecipeValue &right, const Token &token)
  {
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue("binary"), "expression");
    node.push_back(RecipeValue(op), "operator");
    node.push_back(left, "left");
    node.push_back(right, "right");
    add_position(node, token);
    return node;
  }

  RecipeValue additive()
  {
    RecipeValue left = this->multiplicative();
    while (this->is_symbol(this->lexer.peek(), "+") || this->is_symbol(this->lexer.peek(), "-"))
    {
      Token op = this->lexer.next();
      this->skip_newlines_inside();
      left = this->binary(op.text, left, this->multiplicative(), op);
    }
    return left;
  }

  RecipeValue multiplicative()
  {
    RecipeValue left = this->range();
    while (this->is_symbol(this->lexer.peek(), "*") || this->is_symbol(this->lexer.peek(), "/"))
    {
      Token op = this->lexer.next();
      this->skip_newlines_inside();
      left = this->binary(op.text, left, this->range(), op);
    }
    return left;
  }

  RecipeValue range()
  {
    RecipeValue left = this->unary();
    if (this->is_symbol(this->lexer.peek(), ":"))
    {
      Token op = this->lexer.next();
      left = this->binary(":", left, this->unary(), op);
    }
    return left;
  }

  RecipeValue unary()
  {
    const Token &token = this->lexer.peek();
    if (this->is_symbol(token, "-") || this->is_symbol(token, "+"))
    {
      Token op = this->lexer.next();
      RecipeValue operand = this->unary();
      // a negative number literal stays a number
      if (op.text == "-" && operand["expression"].as_string() == "number")
      {
        RecipeValue number = RecipeValue::list();
        number.push_back(RecipeValue("number"), "expression");
        number.push_back(RecipeValue(-operand["value"].as_double()), "value");
        add_position(number, op);
        return number;
      }
      if (op.text == "+")
        return operand;
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("unary"), "expression");
      node.push_back(RecipeValue("-"), "operator");
      node.push_back(operand, "operand");
      add_position(node, op);
      return node;
    }
    return this->power();
  }

  RecipeValue power()
  {
    RecipeValue base = this->postfix();
    if (this->is_symbol(this->lexer.peek(), "^"))
    {
      Token op = this->lexer.next();
      return this->binary("^", base, this->unary(), op);
    }
    return base;
  }

  RecipeValue postfix()
  {
    Token start = this->lexer.peek();
    RecipeValue node = this->primary();
    while (true)
    {
      const Token &token = this->lexer.peek();
      if (this->is_symbol(token, "("))
      {
        if (node["expression"].as_string() != "name")
          this->lexer.error(token.line, token.column, "only a named function can be called.");
        this->lexer.next();
        RecipeValue call = RecipeValue::list();
        call.push_back(RecipeValue("call"), "expression");
        call.push_back(node["name"], "function");
        call.push_back(node["namespace"], "namespace");
        call.push_back(this->arguments(")"), "arguments");
        add_position(call, start);
        node = call;
        continue;
      }
      if (this->is_symbol(token, "$"))
      {
        // a member, as in R: data$y, parameters$theta (in calls to other languages)
        this->lexer.next();
        Token member = this->expect_identifier("a name after \"$\"");
        RecipeValue node_member = RecipeValue::list();
        node_member.push_back(RecipeValue("member"), "expression");
        node_member.push_back(node, "target");
        node_member.push_back(RecipeValue(member.text), "name");
        add_position(node_member, start);
        node = node_member;
        continue;
      }
      if (this->is_symbol(token, "["))
      {
        this->lexer.next();
        RecipeValue indexed = RecipeValue::list();
        indexed.push_back(RecipeValue("index"), "expression");
        indexed.push_back(node, "target");
        indexed.push_back(this->arguments("]"), "indices");
        add_position(indexed, start);
        node = indexed;
        continue;
      }
      break;
    }
    return node;
  }

  // arguments up to the closing symbol; "name = value" is a keyword argument
  RecipeValue arguments(const std::string &close)
  {
    RecipeValue output = RecipeValue::list();
    if (this->is_symbol(this->lexer.peek(), close))
    {
      this->lexer.next();
      return output;
    }
    while (true)
    {
      const Token &token = this->lexer.peek();
      if (token.type == TokenType::Identifier && this->is_symbol(this->lexer.peek(1), "=") && close == ")")
      {
        Token name = this->lexer.next();
        this->lexer.next();
        output.push_back(this->expression(), name.text);
      }
      else
      {
        output.push_back(this->expression());
      }
      Token separator = this->lexer.next();
      if (this->is_symbol(separator, close))
        break;
      if (!this->is_symbol(separator, ","))
        this->lexer.error(separator.line, separator.column,
                          "expected \",\" or \"" + close + "\", found " + describe(separator) + ".");
    }
    return output;
  }

  RecipeValue primary()
  {
    Token token = this->lexer.next();
    RecipeValue node = RecipeValue::list();
    switch (token.type)
    {
      case TokenType::Number:
        node.push_back(RecipeValue("number"), "expression");
        node.push_back(RecipeValue(token.number), "value");
        if (!token.text.empty() && token.text.back() == 'L')
          node.push_back(RecipeValue(true), "integer");
        add_position(node, token);
        return node;
      case TokenType::String:
        node.push_back(RecipeValue("string"), "expression");
        node.push_back(RecipeValue(token.text), "value");
        add_position(node, token);
        return node;
      case TokenType::Identifier:
      {
        const std::string &t = token.text;
        if (t == "TRUE" || t == "true" || t == "FALSE" || t == "false")
        {
          node.push_back(RecipeValue("logical"), "expression");
          node.push_back(RecipeValue(t == "TRUE" || t == "true"), "value");
          add_position(node, token);
          return node;
        }
        if (t == "Inf" || t == "NaN")
        {
          node.push_back(RecipeValue("number"), "expression");
          node.push_back(RecipeValue(t == "Inf" ? std::numeric_limits<double>::infinity()
                                                : std::numeric_limits<double>::quiet_NaN()), "value");
          add_position(node, token);
          return node;
        }
        if (this->is_symbol(this->lexer.peek(), "::"))
        {
          this->lexer.next();
          Token name = this->expect_identifier("a name after \"" + t + "::\"");
          RecipeValue named = this->name_node(name, t);
          // the position of the whole name
          named.set("line", RecipeValue(int(token.line)));
          named.set("column", RecipeValue(int(token.column)));
          return named;
        }
        return this->name_node(token, "");
      }
      case TokenType::Symbol:
        if (token.text == "(")
        {
          RecipeValue inner = this->expression();
          this->expect_symbol(")", "to close \"(\"");
          return inner;
        }
        if (token.text == "[")
        {
          node.push_back(RecipeValue("vector"), "expression");
          node.push_back(this->arguments("]"), "elements");
          add_position(node, token);
          return node;
        }
        break;
      default:
        break;
    }
    this->lexer.error(token.line, token.column, "expected a value, found " + describe(token) + ".");
  }

};

}

RecipeValue parse_ilike_language(const std::string &text, const std::string &file_name)
{
  Parser parser(text, file_name);
  return parser.program();
}

bool is_ilike_language(const std::string &text)
{
  size_t i = 0;
  while (i < text.size())
  {
    if (std::isspace((unsigned char)text[i])) { ++i; continue; }
    if (text.compare(i, 2, "//") == 0)
    {
      while (i < text.size() && text[i] != '\n') ++i;
      continue;
    }
    if (text.compare(i, 2, "/*") == 0)
    {
      size_t end = text.find("*/", i + 2);
      if (end == std::string::npos)
        return false;
      i = end + 2;
      continue;
    }
    break;
  }
  if (text.compare(i, 5, "ilike") != 0)
    return false;
  i += 5;
  if (i >= text.size() || !(text[i] == ' ' || text[i] == '\t'))
    return false;
  while (i < text.size() && (text[i] == ' ' || text[i] == '\t')) ++i;
  return (i < text.size()) && std::isdigit((unsigned char)text[i]);
}

}
