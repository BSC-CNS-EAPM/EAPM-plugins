const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const ReactRefreshWebpackPlugin = require('@pmmmwh/react-refresh-webpack-plugin');

module.exports = (env, argv) => {

  const mode = argv.mode

  return {
    mode: mode,
    entry: {
      docking_analysis: "./src/DockingAnalysis/index.tsx",
    },
    output: {
      path: path.resolve("EAPM", "Pages"),
      filename: "bundle-[name].js",
    },
    devServer: {
      static: "./dist",
      port: 1234,
      historyApiFallback: {
        index: "/docking_analysis.html"
      },
      proxy: [
        {
          context: ['/api'],
          target: 'http://localhost:3000/plugins/pages/eapm.docking_analysis/',
        },
      ],
    },
    infrastructureLogging: {
      level: mode === 'development' ? 'verbose' : 'info',
      debug:
        mode === 'development'
          ? [(name) => name.includes('webpack-dev-server')]
          : false,
    },
    module: {
      rules: [
        {
          test: /\.(js|jsx)$/,
          exclude: /node_modules/,
          use: {
            loader: "babel-loader",
            options: {
              plugins: mode === "development" ? ["react-refresh/babel"] : undefined,
            },
          },
        },
        {
          test: /\.tsx?$/,
          use: [
            {
              loader: "ts-loader",
              options: {
                transpileOnly: true,
              },
            },
          ],
          exclude: /node_modules/,
        },
        {
          test: /\.(sass|less|css|scss)$/,
          use: [
            MiniCssExtractPlugin.loader,
            "css-loader",
            "postcss-loader",
            "sass-loader",
          ],
        },
        {
          test: /\.svg$/i,
          issuer: /\.[jt]sx?$/,
          use: ["@svgr/webpack"],
        },
        {
          test: [/\.bmp$/, /\.gif$/, /\.jpe?g$/, /\.png$/],
          loader: require.resolve("url-loader"),
          options: {
            limit: 10000,
            name: "static/media/[name].[hash].[ext]",
          },
        },
      ],
    },
    resolve: {
      extensions: [".js", ".jsx", ".tsx", ".ts"],
    },
    plugins: [
      new ReactRefreshWebpackPlugin(),
      new MiniCssExtractPlugin(),
      new HtmlWebpackPlugin({
        hash: true,
        title: "Docking Analysis",
        filename: "docking_analysis.html",
      }),
    ],
  }
};

