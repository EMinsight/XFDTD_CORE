#include <iostream>
#include "toml++/toml.hpp"

int main() {
  auto data = toml::parse("config.toml");
  auto title = data["title"].value_or("Untitled");
  auto owner = data["owner"].value_or("nobody");
  auto database = data["database"].value_or("data.db");
  std::cout << "Title: " << title << '\n'
            << "Owner: " << owner << '\n'
            << "Database: " << database << '\n';
}