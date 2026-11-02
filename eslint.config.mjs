import js from "@eslint/js";
import globals from "globals";

export default [
  js.configs.recommended,
  {
    files: ["web/**/*.js"],
    languageOptions: { sourceType: "module", globals: { ...globals.browser, ...globals.worker } },
    rules: {
      "no-unused-vars": ["error", { caughtErrors: "none" }],
      eqeqeq: "error",
      "no-implicit-globals": "error",
    },
  },
];
